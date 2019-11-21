# ppvcf

```
git clone --recursive https://github.com/ldenti/ppvcf.git
cd ppvcf/htslib
make
cd ..
make
export LD_LIBRARY_PATH=./htslib
./main small.vcf.gz 2 0
```
