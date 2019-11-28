# ppvcf
`ppvcf` is a fast c++ library for the parallel parsing of huge vcf files. `ppvcf` is a C++ wrapper for `htslib`: it builds upon `htslib` providing a way to parse in parallel big vcf files containing a huge quantity of samples.

### Installation
`ppvcf` comes as a single `.hpp` header that you can easily include in your source. It depends on:
* [htslib](https://www.htslib.org/)
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
If you have `htslib` installed in your system, to compile the example program, you just have to link against it with the `-lhts` flag:
```bash
g++ -Wall -O3 -std=c++11 example.cpp -o example -lhts -fopenmp
```
Otherwise, if you have installed `htslib` locally, you can compile the example program with:
```bash
g++ -Wall -O3 -std=c++11 -I./htslib example.cpp -o example -L./htslib -lhts -fopenmp
```
Then you can run the example code with:
```
./example tiny.vcf 2
```

### Experiments
To run the experiments:
```
git clone --recursive https://github.com/ldenti/ppvcf.git
cd ppvcf/htslib
make
cd ..
make
export LD_LIBRARY_PATH=./htslib
./main small.vcf.gz 2 0
```