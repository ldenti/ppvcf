# ppvcf

```
git clone --recursive https://github.com/ldenti/ppvcf.git
cd htslib
make
cd ..
g++ -g -O3 -Wall -std=c++11 -I. -I./htslib/htslib main.cpp -o main -L./htslib -lhts -fopenmp
export LD_LIBRARY_PATH=./htslib
./parse_vcf small.vcf.gz 2 0
```
