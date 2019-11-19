# parse_vcf

```
git clone --recursive https://github.com/ldenti/parse_vcf.git
cd htslib
make
cd ..
g++ -g -O0 -Wall -std=c++11 -I. -I./htslib/htslib main.cpp -o main -L./htslib -lhts -fopenmp
export LD_LIBRARY_PATH=./htslib
./parse_vcf small.vcf 2
```
