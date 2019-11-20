#include <iostream>

#include "vcf.h"

#include "vcf_file.hpp"

using namespace std;

int parse_vcf_htslib_ALL(char* vcf_path) {
  htsFile *vcf = bcf_open(vcf_path, "r");
  bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
  bcf1_t *vcf_record = bcf_init();

  int n_individuals = bcf_hdr_nsamples(vcf_header);
  
  while (bcf_read(vcf, vcf_header, vcf_record) == 0) {
    bcf_unpack(vcf_record, BCF_UN_ALL);
    int32_t *gt_arr = NULL, ngt = 0;
    int ngt_ret_value = bcf_get_genotypes(vcf_header, vcf_record, &gt_arr, &ngt);

    if ( ngt<=0 ) return 1; // GT not present

    int ploidy = ngt / n_individuals;
    for(int i=0; i<n_individuals; ++i) {
      int32_t *curr_gt = gt_arr + i * ploidy;
      int all_1;
      int all_2;
      if (curr_gt[1] == bcf_int32_vector_end) {
        all_1 = bcf_gt_allele(curr_gt[0]);
        all_2 = bcf_gt_allele(curr_gt[0]);
      } else {
	cout << curr_gt[0] << " " << bcf_gt_allele(curr_gt[0]) << endl;
        all_1 = bcf_gt_allele(curr_gt[0]);
        all_2 = bcf_gt_allele(curr_gt[1]);
      }
      if(all_1 < 0) all_1 = 0;
      if(all_2 < 0) all_2 = 0;

      GT gt(all_1, all_2, true);
    }
    free(gt_arr);
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

  if(mode == 0) {
    VCF vcf (vcf_path, n_threads);
    while(true) {
      bool a = vcf.parse(100);
      if(!a) {
	break;
      }
    }
  } else
    parse_vcf_htslib_ALL(vcf_path);
  return 0;
}
