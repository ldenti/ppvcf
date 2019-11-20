#ifndef _VARIANT_HPP_
#define _VARIANT_HPP_

#include <numeric>
#include <string>
#include <vector>

using namespace std;

struct GT {
  uint8_t a1, a2;
  bool phased;

  GT(uint8_t a1_ = -1, uint8_t a2_ = -1, bool phased_ = true) {
    a1 = a1_;
    a2 = a2_;
    phased = phased_ || a1 == a2;
  }
  ~GT() { }

  string to_str() const {
    return phased ? to_string(a1)+"|"+to_string(a2) : to_string(a1)+"/"+to_string(a2);
  }
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
  
  Variant(const uint32_t _nsamples) :
    nsamples(_nsamples),
    gti(0),
    genotypes(nsamples) {  }
  ~Variant() {  }

  void update_till_info(bcf_hdr_t *header, bcf1_t *record) {
    seq_name = bcf_hdr_id2name(header, record->rid);
    ref_pos = record->pos;
    idx = record->d.id;
    ref_sub = record->d.allele[0];

    for (int i = 1; i < record->n_allele; ++i) {
      char *curr_alt = record->d.allele[i];
      if (curr_alt[0] != '<')
        alts.push_back(string(curr_alt));
    }

    quality = record->qual;
    filter = "PASS"; // TODO: get filter string
    info = ".";      // TODO: get info string
  }

  void add_genotype(const GT& gt) {
    assert(gti < nsamples);
    genotypes[gti++] = gt;
  }

  /**
   * Deleted functions (see MALVA):
   * - extract_frequencies from INFO
   * - get i-th allele
   **/
};

#endif
