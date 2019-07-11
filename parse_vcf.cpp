#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <omp.h>

#include "vcf.h"
#include "kstring.h"
#include "variant.hpp"

using namespace std;

struct GenoType {
  int a1, a2;
  bool phased;

  GenoType() {
    a1 = 0;
    a2 = 0;
    phased = true;
  }

  GenoType(int a1_, int a2_, bool phased_) {
    a1 = a1_;
    a2 = a2_;
    phased = phased_ || a1 == a2;
    cout << a1 << "-" << a2 << "-" << phased << endl;
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

GenoType extract_genotype(char* gt) {
  cout << gt << endl;
  if(!strcmp(gt, ".")) {
    return GenoType();
  }
  bool phased = strchr(gt, '/') != NULL;
  int all1, all2;
  char *all = strtok (gt,"/|");
  all1 = atoi(all);
  all = strtok(NULL,"/|");
  all2 = all==NULL ? all1 : atoi(all);
  cout << "@ " << all1 << "-" << all2 << "-" << phased << endl;
  return GenoType(all1, all2, phased);
}

void parse_format_field(char *vcf_line) {
  char *samples = get_samples(vcf_line);
  
  for(;;) {
    char *next_tab = strchr(samples, '\t');
    char *gt = strtok(samples,":\t");
    GenoType g = extract_genotype(gt);
    if(next_tab == NULL) break;
    samples = next_tab+1;
  }
  // while (sample != NULL) {
  //   cout << strlen(sample) << " " << sample << endl;
    
  //   // char *gt = strtok(sample,":");
  //   // cout << gt << endl;

  //   sample = strtok(NULL,"\t");
  // }
}

void analyze(vector<kstring_t*> &vcf_lines, vector<bcf1_t> &vcf_records, bcf_hdr_t *vcf_header, const int &n_threads) {
#pragma omp parallel for num_threads (n_threads) shared (vcf_header, vcf_lines, vcf_records)
  for(int i=0; i<(int)vcf_lines.size(); ++i) {
    parse_format_field(vcf_lines[i]->s);
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
