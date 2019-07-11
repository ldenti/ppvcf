#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <omp.h>

#include "vcf.h"
#include "kstring.h"
// #include "variant.hpp"

using namespace std;

struct GT { // This name is already used in variant.hpp
  uint8_t a1, a2;
  bool phased;

  GT() {
    a1 = 0;
    a2 = 0;
    phased = true;
  }

  GT(uint8_t a1_, uint8_t a2_, bool phased_) {
    a1 = a1_;
    a2 = a2_;
    phased = phased_ || a1 == a2;
  }

  string to_str() {
    return phased ? to_string(a1)+"|"+to_string(a2) : to_string(a1)+"/"+to_string(a2);
  }

  void print() {
    if(phased)
      cout << (int)a1 << "|" << (int)a2 << endl;
    else
      cout << (int)a1 << "/" << (int)a2 << endl;
  }
};

char* get_samples(char *vcf_line) {
  int offset = 0;
  char* fmt_suffix;
  int fmt_size = 0;
  bool fmt_flag = false;
  while(!fmt_flag) {
    offset += strlen(vcf_line + offset) + 1;
    /**
     * From VCF format specification: "The first sub-field must always
     * be the genotype (GT) if it is present. There are no required
     * sub-fields."
     **/
    if(strncmp(vcf_line + offset, "GT\t", 3) == 0 || strncmp(vcf_line + offset, "GT:", 3) == 0) {
      fmt_flag = true;
      fmt_suffix = vcf_line + offset;
      char* fmt_end = strchr(fmt_suffix, '\t');
      if(fmt_end==NULL) {
	cerr << "Error - No samples" << endl;
	exit(1);
      }
      fmt_size = fmt_end - fmt_suffix;
    }
  }
  return vcf_line + offset + fmt_size + 1;
}

GT extract_genotype(char* gt) {
  if(!strcmp(gt, ".")) {
    return GT();
  }
  bool phased = strchr(gt, '/') != NULL;
  int all1, all2;
  char *all = strtok (gt,"/|");
  all1 = atoi(all);
  all = strtok(NULL,"/|");
  all2 = all==NULL ? all1 : atoi(all);
  return GT(all1, all2, phased);
}

void parse_format_field(char *vcf_line, const int &n_samples) {
  char *samples = get_samples(vcf_line);

  vector<GT> genotypes (n_samples);
  int i = 0;
  for(;;) {
    char *next_tab = strchr(samples, '\t');
    char *gtc = strtok(samples,":\t");
    GT gt = extract_genotype(gtc);
    genotypes[i] = gt;
    if(next_tab == NULL) break;
    samples = next_tab+1;
  }
}

void analyze(vector<kstring_t*> &vcf_lines, const int &n_samples, const int &n_threads) {
#pragma omp parallel for num_threads (n_threads) shared (vcf_lines, n_samples)
  for(int i=0; i<(int)vcf_lines.size(); ++i) {
    parse_format_field(vcf_lines[i]->s, n_samples);
    // Variant v(vcf_header, &vcf_records[i], "AF");
    // cerr << v.idx << " " << v.ref_pos << endl;
  }
}

int parse_vcf_htslib_mixed(char* vcf_path, const int &n_threads) {
  htsFile *vcf = bcf_open(vcf_path, "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  int n_samples = bcf_hdr_nsamples(vcf_header);
  bcf1_t *vcf_record = bcf_init();
  vcf_record->max_unpack = BCF_UN_INFO;

  vector<kstring_t*> vcf_lines;
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_INFO); // here we unpack till INFO field

    // Maybe here we can just copy vcf->line.s and make a vector of char*?
    kstring_t *ks = new kstring_t();
    ks_initialize(ks);
    kputsn(vcf->line.s, vcf->line.l, ks);
    vcf_lines.push_back(ks);

    if(vcf_lines.size() == 20000) {
      analyze(vcf_lines, n_samples, n_threads);
      for(kstring_t* ks : vcf_lines)
	ks_free(ks);
      vcf_lines.clear();
    }
  }
  if(!vcf_lines.empty()) {
    analyze(vcf_lines, n_samples, n_threads);
    vcf_lines.empty();
  }

  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  return 0;
}

int parse_vcf_htslib_ALL(char* vcf_path, const int &n_threads) {
  htsFile *vcf = bcf_open(vcf_path, "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  bcf1_t *vcf_record = bcf_init();
  vcf_record->max_unpack = BCF_UN_ALL;

  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_ALL);
  }

  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  return 0;
}

int main(int argc, char *argv[]) {
  char* vcf_path = argv[1];
  int n_threads = atoi(argv[2]);
  int mode = atoi(argv[3]);
  if(mode == 0)
    parse_vcf_htslib_mixed(vcf_path, n_threads);
  else
    parse_vcf_htslib_ALL(vcf_path, n_threads);
  // cerr << "END" << endl;
  return 0;
}
