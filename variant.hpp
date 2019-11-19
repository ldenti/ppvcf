#ifndef _VARIANT_HPP_
#define _VARIANT_HPP_

#include <numeric>
#include <string>
#include <vector>

using namespace std;

struct GT {
  uint8_t a1, a2;
  bool phased;

  GT(uint8_t a1_ = 0, uint8_t a2_ = 0, bool phased_ = true) {
    a1 = a1_;
    a2 = a2_;
    phased = phased_ || a1 == a2;
  }

  string to_str() {
    return phased ? to_string(a1)+"|"+to_string(a2) : to_string(a1)+"/"+to_string(a2);
  }
};

class Variant {
private:
  string seq_name;
  int ref_pos;                                // Variant position 0-based
  string idx;                            // ID
  string ref_sub;                        // Reference base{s}
  vector<string> alts;              // List of alternatives
  float quality;                              // Quality field
  string filter;                         // Filter field
  string info;                           // Info field
  
  vector<GT> genotypes; // full list of genotypes

public:
  Variant(const uint32_t nsamples) : genotypes(nsamples) { }
  ~Variant() {
  }

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
    filter = "PASS"; // TODO: get filter string from VCF
    info = ".";      // TODO: get info string from VCF
  }

  /**
   * Deleted functions (see MALVA):
   * - extract_frequencies from INFO
   * - get i-th allele
   **/
};

#endif
