#ifndef _VCFFILE_HPP_
#define _VCFFILE_HPP_

#include <omp.h>
#include <tuple>
#include <vector>

#include "htslib/vcf.h"

#include "variant.hpp"

using namespace std;

typedef struct {
  size_t size;
  char *seq;
} str_w_l;

class VCF {
private:
  htsFile *vcf;
  bcf_hdr_t *vcf_header;
  bcf1_t *vcf_record;

  int n_samples;
  vector<Variant> variants;
  vector<str_w_l> fmt_lines;
  size_t block_size;
  size_t to_parse;

  int n_threads;

public:
  VCF(char *vcf_path, const int nths, const int block_size_);
  ~VCF();

  Variant front();
  Variant back();
  vector<Variant>::iterator begin();
  vector<Variant>::iterator end();
  vector<Variant>::const_iterator begin() const;
  vector<Variant>::const_iterator end() const;
  vector<Variant>::const_iterator cbegin() const;
  vector<Variant>::const_iterator cend() const;

  /**
   * Parse n lines from VCF. Returns 0 if it reads less than
   * n lines (due to EOF).
   **/
  bool parse();

private:
  /**
   * This function cuts the FORMAT description field
   **/
  char *get_samples(char *const fmt_line, const size_t len);

  /**
   * Transform a char* representing a genotype to the struct GT.
   **/
  tuple<uint8_t, uint8_t, bool> extract_genotype(char *start);

  void fill_variant(const uint32_t var_idx);

  void fill_genotypes();
};

#endif
