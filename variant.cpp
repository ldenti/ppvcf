#include "variant.hpp"

/*---------- GT methods ----------*/

GT::GT(uint8_t a1_, uint8_t a2_, bool phased_) : a1(a1_), a2(a2_) {
  phased = phased_ || a1 == a2;
}

GT::~GT() {}

string GT::to_str() const {
  return phased ? to_string(a1) + "|" + to_string(a2)
                : to_string(a1) + "/" + to_string(a2);
}

/*---------- Variant methods ----------*/
Variant::Variant(const uint32_t _nsamples)
    : nsamples(_nsamples), gti(0), genotypes(nsamples) {}

Variant::~Variant() {}

void Variant::update_till_info(bcf_hdr_t *header, bcf1_t *record) {
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

void Variant::add_genotype(const GT &gt) {
  assert(gti < nsamples);
  genotypes[gti++] = gt;
}
