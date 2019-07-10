#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <omp.h>

#include "vcf.h"
#include "kstring.h"
#include "variant.hpp"

using namespace std;

void parse_format_field(kstring_t *vcf_line, const bcf_hdr_t *vcf_header, bcf1_t *vcf_record) {
  int offset = 0;
  bool fmt_flag = false;
  int nfield = 0;
  while(!fmt_flag) {
    ++nfield;
    // cout << vcf_line.s + offset << " ";
    offset += strlen(vcf_line->s + offset) + 1;
    // Here I'm assuming FORMAT field to start with GT
    if(strncmp(vcf_line->s + offset, "GT\t", 3) == 0 || strncmp(vcf_line->s + offset, "GT:", 3) == 0) {
      fmt_flag = true;
      char* fmt_suffix = vcf_line->s + offset;
      char* fmt_end = strchr(fmt_suffix, '\t');
      // cerr << "!!! " << fmt_suffix << endl;
      if(fmt_end==NULL) {
	cerr << "!!! " << vcf_line->s + offset << endl;
	cerr << "Error - Invalid GT" << endl;
	exit(1);
      }
      int fmt_size = fmt_end - fmt_suffix;
      char *q = fmt_suffix + fmt_size; //end of FORMAT field
      *q = 0;
      vcf_parse_format(vcf_line, vcf_header, vcf_record, fmt_suffix, q);
    }
  }
}

void analyze(vector<kstring_t*> &vcf_lines, vector<bcf1_t> &vcf_records, bcf_hdr_t *vcf_header, const int &n_threads) {
#pragma omp parallel for num_threads (n_threads) shared (vcf_header, vcf_lines, vcf_records)
  for(int i=0; i<(int)vcf_lines.size(); ++i) {
    parse_format_field(vcf_lines[i], vcf_header, &vcf_records[i]);
    // Variant v(vcf_header, &vcf_records[i], "AF");
    // cerr << v.idx << " " << v.ref_pos << endl;
  }
  for(auto &s : vcf_lines) {
    ks_free(s);
  }  
}

int parse_vcf_htslib(char* vcf_path, const int &n_threads) {
  htsFile *vcf = bcf_open(vcf_path, "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  bcf1_t *vcf_record = bcf_init();
  vcf_record->max_unpack = BCF_UN_INFO;

  vector<kstring_t*> vcf_lines;
  vector<bcf1_t> vcf_records;
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_INFO); // here we unpack till INFO field
    vcf_records.push_back(*vcf_record);
    kstring_t *ks = new kstring_t();
    ks_initialize(ks);
    kputsn(vcf->line.s, vcf->line.l, ks);
    vcf_lines.push_back(ks);

    if(vcf_lines.size() == 20000) {
      analyze(vcf_lines, vcf_records, vcf_header, n_threads);
      vcf_lines.clear();
      vcf_records.clear();
    }
  }
  if(!vcf_lines.empty()) {
    analyze(vcf_lines, vcf_records, vcf_header, n_threads);
    vcf_lines.empty();
    vcf_records.empty();
  }

  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);

  return 0;
}

int main(int argc, char *argv[]) {
  char* vcf_path = argv[1];
  int n_threads = atoi(argv[2]);
  parse_vcf_htslib(vcf_path, n_threads);
  cerr << "END" << endl;
  return 0;
}
