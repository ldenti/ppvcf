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
  Variant front() { return variants.front(); }
  Variant back() { return variants.back(); }
  vector<Variant>::iterator begin() { return variants.begin(); }
  vector<Variant>::iterator end() { return variants.end(); }
  vector<Variant>::const_iterator begin() const { return variants.begin(); }
  vector<Variant>::const_iterator end() const { return variants.end(); }
  vector<Variant>::const_iterator cbegin() const { return variants.cbegin(); }
  vector<Variant>::const_iterator cend() const { return variants.cend(); }
  
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
    variants.clear();
    uint32_t i = 0;
    while (i<n && bcf_read(vcf, vcf_header, vcf_record) == 0) {
      // 1. we unpack till INFO field
      bcf_unpack(vcf_record, BCF_UN_INFO);
      variants.push_back(Variant(n_samples));
      variants.back().update_till_info(vcf_header, vcf_record);
      
      // 2. we store the line
      // FIXME: we should store only the FORMAT fields
      char* fmt_line = new char[vcf->line.l+1];
      memcpy(fmt_line, vcf->line.s, vcf->line.l+1);
      fmt_lines.push_back(fmt_line);

      ++i;
    }

    fill_genotypes();

    for(uint i=0; i<fmt_lines.size(); ++i)
      free(fmt_lines[i]);
    fmt_lines.clear();

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
  GT extract_genotype(char* start) {
    if(!strcmp(start, "."))
      return GT();

    int all1, all2;
    char* start2 = start;;

    for(; *start2 != '|' && *start2 != '/'; ++start2) ;
    bool phased = *start2 == '|';
    *start2 = '\0';
    ++start2;

    all1 = atoi(start);
    all2 = atoi(start2);

    // assert(all1 <= alts.size() && all2 <= alts.size());

    return GT(all1, all2, phased);
  }

  void fill_variant(const uint32_t var_idx) {
    char *samples = get_samples(fmt_lines[var_idx]);

    int fmt_field_idx = 0; // if we find a :, then we increment this
    char* start = samples;
    int n = strlen(samples);
    for(int i=0; i<=n; ++i) {
      char c = *(samples+i);
      if(c == '\t' || c == '\0') {
	*(samples+i) = '\0';
	GT gt = extract_genotype(start);
	// cout << gt.to_str() << endl;
	variants[var_idx].add_genotype(gt);
	start = samples+i+1;
      } else if(c == ':') {
	++fmt_field_idx;
	*(samples+i) = '\0';
      } else {  }
    }
  }

  void fill_genotypes() {
#pragma omp parallel for num_threads (n_threads) shared (fmt_lines, variants)
    for(uint i=0; i<fmt_lines.size(); ++i)
      fill_variant(i);
  }
};

#endif
