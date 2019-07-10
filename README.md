# parse_vcf

```
git clone --recursive https://github.com/ldenti/parse_vcf.git
cd htslib
git checkout FORMAT_parse
make
cd ..
g++ -O0 -g -Wall -std=c++11 -I. -I./htslib/htslib parse_vcf.cpp -o parse_vcf -L./htslib -lhts -fopenmp -pthread
export LD_LIBRARY_PATH=./htslib
./parse_vcf small.vcf 2
```
