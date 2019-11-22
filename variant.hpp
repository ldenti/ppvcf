#ifndef _VARIANT_HPP_
#define _VARIANT_HPP_

#include <numeric>
#include <string>
#include <vector>

#include "htslib/vcf.h"

using namespace std;

struct GT {
  uint8_t a1, a2;
  bool phased;

  GT(uint8_t a1_ = -1, uint8_t a2_ = -1, bool phased_ = true);

  ~GT();

  string to_str() const;
};

class Variant {
private:
  string seq_name;
  int ref_pos;
  string idx;
  string ref_sub;
  vector<string> alts;
  float quality;
  string filter;
  string info;

  uint32_t nsamples;
  uint32_t gti;

public:
  vector<GT> genotypes; // full list of genotypes
  
  Variant(const uint32_t _nsamples);
  ~Variant();

  void update_till_info(bcf_hdr_t *header, bcf1_t *record);

  void add_genotype(const GT& gt);

  /**
   * Deleted functions (see MALVA):
   * - extract_frequencies from INFO
   * - get i-th allele
   **/
};

#endif
