#ifndef _VCFFILE_HPP_
#define _VCFFILE_HPP_

#include <vector>
#include <omp.h>

#include "vcf.h"

#include "variant.hpp"

using namespace std;

class VCF {
private:
  htsFile *vcf;
  bcf_hdr_t *vcf_header;
  bcf1_t *vcf_record;

  int n_samples;
  vector<Variant> variants;
  vector<char*> fmt_lines;

  int n_threads;
  
public:
  VCF(char* vcf_path, const int &nths) {
    vcf = bcf_open(vcf_path, "r");
    vcf_header = bcf_hdr_read(vcf);
    vcf_record = bcf_init();
    vcf_record->max_unpack = BCF_UN_INFO;

    n_samples = bcf_hdr_nsamples(vcf_header);

    n_threads = nths;
  }

  ~VCF() {
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf);
  }

  /**
   * Parse n lines from VCF. Returns 0 if it reads less than
   * n lines (due to EOF).
   **/
  bool parse(const uint32_t n) {
    uint32_t i = 0;
    while (i<n && bcf_read(vcf, vcf_header, vcf_record) == 0) {
      // 1. we unpack till INFO field
      bcf_unpack(vcf_record, BCF_UN_INFO);
      variants.push_back(Variant(n_samples));
      variants.back().update_till_info(vcf_header, vcf_record);
      
      // 2. we store the rest of the line (ie FORMAT fields)
      char* fmt_line = new char[vcf->line.l];
      memcpy(fmt_line, vcf->line.s, vcf->line.l);
      fmt_lines.push_back(fmt_line);

      ++i;
    }

    fill_genotypes();

    return i==n;
  }

private:
  /**
   * This function cuts the FORMAT description field
   **/
  char* get_samples(char *fmt_line) {
    int offset = 0;
    char* fmt_suffix;
    int fmt_size = 0;
    bool fmt_flag = false;
    while(!fmt_flag) {
      offset += strlen(fmt_line + offset) + 1;
      /**
       * From VCF format specification: "The first sub-field must always
       * be the genotype (GT) if it is present. There are no required
       * sub-fields."
       **/
      if(strncmp(fmt_line + offset, "GT\t", 3) == 0 || strncmp(fmt_line + offset, "GT:", 3) == 0) {
  	fmt_flag = true;
  	fmt_suffix = fmt_line + offset;
  	char* fmt_end = strchr(fmt_suffix, '\t');
  	if(fmt_end==NULL) {
  	  cerr << "Error - No samples" << endl;
  	  exit(1);
  	}
  	fmt_size = fmt_end - fmt_suffix;
      }
    }
    return fmt_line + offset + fmt_size + 1;
  }


  /**
   * Transform a char* representing a genotype to the struct GT.
   **/
  GT extract_genotype(char* gt) {
    if(!strcmp(gt, ".") or !strcmp(gt, "./.") or !strcmp(gt, ".|."))
      return GT();

    bool phased = strchr(gt, '|') != NULL;
    int all1, all2;
    char *gt_ptr;
    char *all = strtok_r(gt,"/|", &gt_ptr);
    all1 = atoi(all);
    all = strtok_r(NULL,"/|", &gt_ptr);
    all2 = all==NULL ? all1 : atoi(all);

    // assert(all1 <= alts.size() && all2 <= alts.size());

    return GT(all1, all2, phased);
  }

  void fill_variant(const uint32_t i) {
    char *samples = get_samples(fmt_lines[i]);

    int j = 0;
    for(;;) {
      char *next_tab = strchr(samples, '\t');
      char *gt_ptr;
      char *gtc = strtok_r(samples,":\t", &gt_ptr);
      GT gt = extract_genotype(gtc);
      // variants[i].add_genotype(gt);
      if(next_tab == NULL) break;
      samples = next_tab+1;
      ++j;
    }
  }
  
  void fill_genotypes() {
#pragma omp parallel for num_threads (n_threads) shared (fmt_lines, variants)
    for(int i=0; i<(int)fmt_lines.size(); ++i) {
      fill_variant(i);
      free(fmt_lines[i]);
    }
    fmt_lines.clear();
  }
};

#endif
