#include <iostream>
#include "ppvcf.hpp"

int main(int argc, char *argv[]) {
  char *vcf_path = argv[1];
  int n_threads = atoi(argv[2]);

  // Initialize a VCF object
  VCF vcf(vcf_path, n_threads, 10000);

  // Read the file in blocks of 10000 variants
  while (vcf.parse()) {
    // Do something with the variants
    for(const Variant v : vcf.get_variants())
      std::cout << v.get_pos() << std::endl;
  }
  // Manage the last block of variants
  for(const Variant v : vcf.get_variants())
    std::cout << v.get_pos() << std::endl;

  return 0;
}
