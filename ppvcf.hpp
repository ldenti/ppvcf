#ifndef _PPVCF_HPP_
#define _PPVCF_HPP_

#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>
#include <map>

#include <omp.h>

#include "htslib/vcf.h"

using namespace std;

/*---------- Auxiliary structs ----------*/
typedef struct {
  size_t size;
  char *seq;
} str_w_l;

/*---------- Genotype struct ----------*/
struct GT {
  uint8_t a1, a2;
  bool phased;

  GT(uint8_t a1_ = -1, uint8_t a2_ = -1, bool phased_ = true) : a1(a1_), a2(a2_) {
    phased = phased_ || a1 == a2;
  }
  ~GT() {}
};

/*---------- Variant class ----------*/
class Variant {
private:
  string seq_name;
  int ref_pos;
  string idx;
  string ref_sub;
  vector<string> alts;
  string quality;
  string filter;
  map<string, string> info;

  uint32_t nsamples;
  uint32_t gti;
  vector<GT> genotypes; // full list of genotypes

  void store_filter(bcf_hdr_t *header, bcf1_t *record) {
    // TODO: store filter differently (not as a single string) if we want to allow filtering
    filter = "";
    if (record->d.n_flt) {
      for (int i = 0; i < record->d.n_flt; ++i) {
	if(i)
	  filter+=";";
	filter += header->id[BCF_DT_ID][record->d.flt[i]].key;
      }
    } else
      filter = ".";
  }

  void store_info(bcf_hdr_t *header, bcf1_t *record) {
    for (int i = 0; i < record->n_info; ++i) {
      bcf_info_t *info_field = &record->d.info[i];
      int key_idx = info_field->key;
      int type = info_field->type;
      const char* key = header->id[BCF_DT_ID][key_idx].key;

      int ndst = 0;
      string value = "";
      if (type == BCF_BT_NULL) {
	bool v = bcf_get_info_flag(header, record, key, &v, &ndst);
	value = to_string(v);
	// FIXME: remember to manage booleans (eg when producing output)
      } else if (type == BCF_BT_INT8 || type == BCF_BT_INT16 || type == BCF_BT_INT32) {
	int *v = NULL;
	bcf_get_info_int32(header, record, key, &v, &ndst);
	value = to_string(*v);
	free(v);
      } else if (type == BCF_BT_FLOAT) {
	float *v = NULL;
	bcf_get_info_float(header, record, key, &v, &ndst);
	value = to_string(*v);
	free(v);
      } else if (type == BCF_BT_CHAR) {
	char *v = NULL;
	bcf_get_info_string(header, record, key, &v, &ndst);
	value = string(v);
	free(v);
      } else {
	cerr << "Unknown type " << type << " (field " << key << ")" << endl;
	exit(1);
      }

      info[key] = value;

      /** // Maybe we can adapt this
	  #define BRANCH(type_t, bcf_ht_t) {				\
	  type_t *value = NULL;						\
	  bcf_get_info_values(header, record, key, (void**)(&value), &ndst, bcf_ht_t); \
	  cout << *value << endl;					\
	  }
	  switch(type) {
	  case BCF_BT_INT8: BRANCH(int, BCF_HT_INT); break;
	  case BCF_BT_INT16: BRANCH(int, BCF_HT_INT); break;
	  case BCF_BT_INT32: BRANCH(int, BCF_HT_INT); break;
	  case BCF_BT_FLOAT: BRANCH(float, BCF_HT_REAL); break;
	  case BCF_BT_CHAR: BRANCH(string, BCF_HT_STR); break;
	  default: cerr << "Unknown type " << type << endl; exit(1);
	  }
	  #undef BRANCH
      **/
    }
  }

public:
  Variant(const uint32_t _nsamples)
    : nsamples(_nsamples), gti(0), genotypes(_nsamples) {}
  ~Variant() {}

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

    quality = to_string(record->qual);
    quality = quality=="nan" ? "." : quality;
    store_filter(header, record);
    store_info(header, record);
  }

  void add_genotype(const uint8_t a1, const uint8_t a2,
		    const bool phased) {
    assert(gti < nsamples);
    genotypes[gti].a1 = a1;
    genotypes[gti].a2 = a2;
    genotypes[gti].phased = phased;
    ++gti;
  }

  string get_info(const string &key) const {
    return info.at(key);
  }
};

/*---------- VCF class ----------*/
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
  VCF(char *vcf_path, const int nths, const int block_size_) {
    vcf = bcf_open(vcf_path, "r");
    vcf_header = bcf_hdr_read(vcf);
    vcf_record = bcf_init();
    vcf_record->max_unpack = BCF_UN_INFO;

    n_samples = bcf_hdr_nsamples(vcf_header);

    n_threads = nths;

    block_size = block_size_;
    to_parse = 0;

    for (size_t i = 0; i < block_size; ++i) {
      str_w_l s;
      s.size = 1024 * 10;
      s.seq = (char *)calloc(s.size, sizeof(char));

      fmt_lines.push_back(s);
    }
  }
  ~VCF() {
    for (size_t i = 0; i < block_size; ++i) {
      free(fmt_lines[i].seq);
    }
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf);
  }

  /**
   * Parse n lines from VCF. Returns 0 if it reads less than
   * n lines (due to EOF).
   **/
  bool parse() {
    variants.clear();
    uint32_t i = 0;
    uint32_t n = block_size;
    while (i < block_size && bcf_read(vcf, vcf_header, vcf_record) == 0) {
      // 1. we unpack till INFO field
      bcf_unpack(vcf_record, BCF_UN_INFO);
      variants.push_back(Variant(n_samples));
      variants.back().update_till_info(vcf_header, vcf_record);

      // 2. we store the line
      char *samples = get_samples(vcf->line.s, vcf->line.l + 1);
      size_t samples_s = vcf->line.l - (samples - vcf->line.s);
      while (samples_s > fmt_lines[i].size) {
	fmt_lines[i].size *= 1.6;
	fmt_lines[i].seq = (char *)realloc(fmt_lines[i].seq, fmt_lines[i].size);
      }
      strncpy(fmt_lines[i].seq, samples, samples_s);
      ++i;
    }
    to_parse = i;
    fill_genotypes();

    return i == n;
  }

  vector<Variant> get_variants() const {
    return variants;
  }

private:
  /**
   * This function cuts the FORMAT description field
   **/
  char *get_samples(char *const fmt_line, const size_t len) {
    /**
     * From VCF format specification: "The first sub-field must always
     * be the genotype (GT) if it is present. There are no required
     * sub-fields."
     **/
    char *gt_start = fmt_line;
    char *const gt_end = fmt_line + len;
    for (; gt_start != gt_end && (*gt_start != 'G' || *(gt_start + 1) != 'T' ||
				  *(gt_start + 2) != '\t');
	 ++gt_start)
      ;
    return gt_start;
  }

  /**
   * Transform a char* representing a genotype to the struct GT.
   **/
  tuple<uint8_t, uint8_t, bool> extract_genotype(char *start) {
    if (*start == '.')
      return tuple<uint8_t, uint8_t, bool>(-1, -1, false);

    uint8_t all1, all2;
    char *start2 = start;

    for (; *start2 != '|' && *start2 != '/'; ++start2)
      ;
    bool phased = *start2 == '|';
    *start2 = '\0';

    all1 = atoi(start);
    all2 = atoi(start2 + 1);

    *start2 = '|'; // don't care, just fix

    return std::move(tuple<uint8_t, uint8_t, bool>(all1, all2, phased));
  }

  void fill_variant(const uint32_t var_idx) {
    char *samples = fmt_lines[var_idx].seq;

    // Drop "GT" at the beginning of the line
    for (; *samples != '\t'; ++samples)
      ;
    ++samples;

    char *start = samples;
    char *end = start;
    bool last = *end == '\0';
    tuple<uint8_t, uint8_t, bool> gt;
    while (!last) {
      for (; *end != '\t' && *end != '\0'; ++end)
	;
      last = *end == '\0';
      *end = '\0';
      gt = extract_genotype(start);
      variants[var_idx].add_genotype(get<0>(gt), get<1>(gt), get<2>(gt));
      *end = last ? '\0' : '\t';
      start = end + 1;
      end = start;
    }
    // TODO: manage ":"
  }

  void fill_genotypes() {
#pragma omp parallel for num_threads(n_threads) shared(fmt_lines, variants)
    for (uint i = 0; i < to_parse; ++i) {
      fill_variant(i);
    }
  }

};

#endif
