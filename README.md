# ppvcf
`ppvcf` is a fast C++ library for the parallel parsing of huge vcf files. It provides a way to parse in parallel big VCF files containing a huge quantity of samples.

### Installation
`ppvcf` comes as a single `.hpp` header that you can easily include in your source. It depends on:
* [openmp](https://www.openmp.org/)


### Small example
See `example.cpp` code:
```cpp
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
```
To compile the example program:
```bash
g++ -Wall -O3 -std=c++11 example.cpp -o example -fopenmp
```
To run the example code with 2 threads:
```
./example tiny.vcf 2
```
