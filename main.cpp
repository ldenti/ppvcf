#include <iostream>

#include "vcf_file.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  char* vcf_path = argv[1];
  int n_threads = atoi(argv[2]);

  VCF vcf (vcf_path, n_threads);
  while(true) {
    bool a = vcf.parse(100);
    if(!a) {
      break;
    }
  }
  return 0;
}
