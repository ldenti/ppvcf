#ifndef _PPVCF_HPP_
#define _PPVCF_HPP_

#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#include <omp.h>

#include "htslib/vcf.h"

/*---------- Auxiliary structs ----------*/
//!  Auxiliary struct
/*!
  Auxiliary struct used to store a cstring and its size.
*/
typedef struct {
  /// cstring size
  size_t size;
  /// cstring
  char *seq;
} str_w_l;

/*---------- Genotype struct ----------*/
//! Struct for a single genotype
/*!
  Genotype information of a single sample associated to a variant.
*/
struct GT {
  /// First allele
  uint8_t a1;
  /// Second allele
  uint8_t a2;
  /// Phased flag
  bool phased;

  /**
   * Constructor
   *
   * @param _a1 is the first allele of the genotype
   * @param _a2 is the second allele of the genotype
   * @param _phased is true if the genotype is phased, false otherwise
   **/
  GT(uint8_t _a1 = 0, uint8_t _a2 = 0, bool _phased = true) : a1(_a1), a2(_a2) {
    phased = _phased || a1 == a2;
  }
  /**
   * Destructor
   **/
  ~GT() {}
};

/*---------- Variant class ----------*/
//! Class for a single variant
/*!
  Single line extracted from a VCF file.
*/
class Variant {
private:
  /// Chromosome
  std::string seq_name;
  /// Position (0-based)
  int pos;
  /// Identifier
  std::string idx;
  /// Reference base(s)
  std::string ref;
  /// Alternate base(s)
  std::vector<std::string> alts;
  /// Phred-scaled quality score
  std::string quality;
  /// Filter status
  std::string filter;
  /// Additional information
  std::map<std::string, std::string> info;

  /// Total number of samples
  uint32_t nsamples;
  /// Genotype counter represeting how many GTs have been stored
  uint32_t gti;
  /// Genotype fields
  std::vector<GT> genotypes;

  /**
   * Extract the filter field from record and store it as a string.
   *
   * @param header is the htslib header
   * @param record is a htslib record
   **/
  void store_filter(bcf_hdr_t *header, bcf1_t *record) {
    // TODO: store filter differently (not as a single string) if we want to
    // allow filtering
    filter = "";
    if (record->d.n_flt) {
      for (int i = 0; i < record->d.n_flt; ++i) {
        if (i)
          filter += ";";
        filter += header->id[BCF_DT_ID][record->d.flt[i]].key;
      }
    } else
      filter = ".";
  }

  /**
   * Extract the additional information from record and store them as a
   *map<string, string>.
   *
   * @param header is the htslib header
   * @param record is a htslib record
   **/
  void store_info(bcf_hdr_t *header, bcf1_t *record) {
    for (int i = 0; i < record->n_info; ++i) {
      bcf_info_t *info_field = &record->d.info[i];
      int key_idx = info_field->key;
      int type = info_field->type;
      const char *key = header->id[BCF_DT_ID][key_idx].key;

      int ndst = 0;
      std::string value = "";
      if (type == BCF_BT_NULL) {
        bool v = bcf_get_info_flag(header, record, key, &v, &ndst);
        value = std::to_string(v);
        // FIXME: remember to manage booleans (eg when producing output)
      } else if (type == BCF_BT_INT8 || type == BCF_BT_INT16 ||
                 type == BCF_BT_INT32) {
        int *v = NULL;
        bcf_get_info_int32(header, record, key, &v, &ndst);
        value = std::to_string(*v);
        free(v);
      } else if (type == BCF_BT_FLOAT) {
        float *v = NULL;
        bcf_get_info_float(header, record, key, &v, &ndst);
        value = std::to_string(*v);
        free(v);
      } else if (type == BCF_BT_CHAR) {
        char *v = NULL;
        bcf_get_info_string(header, record, key, &v, &ndst);
        value = std::string(v);
        free(v);
      } else {
        cerr << "Unknown type " << type << " (field " << key << ")" << endl;
        exit(1);
      }

      info[key] = value;

      /* // Maybe we can adapt this
         #define BRANCH(type_t, bcf_ht_t) {				\
         type_t *value = NULL;						\
         bcf_get_info_values(header, record, key, (void**)(&value), &ndst,
         bcf_ht_t); \
         cout << *value << endl;					\
         }
         switch(type) {
         case BCF_BT_INT8: BRANCH(int, BCF_HT_INT); break;
         case BCF_BT_INT16: BRANCH(int, BCF_HT_INT); break;
         case BCF_BT_INT32: BRANCH(int, BCF_HT_INT); break;
         case BCF_BT_FLOAT: BRANCH(float, BCF_HT_REAL); break;
         case BCF_BT_CHAR: BRANCH(std::string, BCF_HT_STR); break;
         default: cerr << "Unknown type " << type << endl; exit(1);
         }
         #undef BRANCH
      */
    }
  }

public:
  /**
   * Constructor
   *
   * @param _nsamples is the number of samples
   **/
  Variant(const uint32_t _nsamples)
      : nsamples(_nsamples), gti(0), genotypes(_nsamples) {}
  /**
   * Destructor
   **/
  ~Variant() {}

  /**
   * Extract and store the first 8 fields of record.
   *
   * @param header is the htslib header
   * @param record is a htslib record
   **/
  void update_till_info(bcf_hdr_t *header, bcf1_t *record) {
    seq_name = bcf_hdr_id2name(header, record->rid);
    pos = record->pos;
    idx = record->d.id;
    ref = record->d.allele[0];

    for (int i = 1; i < record->n_allele; ++i) {
      char *curr_alt = record->d.allele[i];
      if (curr_alt[0] != '<')
        alts.push_back(std::string(curr_alt));
    }

    quality = std::to_string(record->qual);
    quality = quality == "nan" ? "." : quality;
    store_filter(header, record);
    store_info(header, record);
  }

  /**
   * Store a single genotype information.
   *
   * @param a1 is the first allele
   * @param a2 is the second allele
   * @param phased is true if the genotype is phased, false otherwise
   **/
  void add_genotype(const uint8_t a1, const uint8_t a2, const bool phased) {
    assert(gti < nsamples);
    genotypes[gti].a1 = a1;
    genotypes[gti].a2 = a2;
    genotypes[gti].phased = phased;
    ++gti;
  }

  /**
   * @return the chromosome of the variant
   **/
  std::string get_seqname() const { return seq_name; }

  /**
   * @return the (0-based) position of the variant
   **/
  int get_pos() const { return pos; }

  /**
   * @return the variant identifier
   **/
  std::string get_idx() const { return idx; }

  /**
   * @return the reference allele
   **/
  std::string get_ref() const { return ref; }

  /**
   * @return the alternate alleles
   **/
  std::vector<std::string> get_alts() const { return alts; }

  /**
   * @return the i-th (1-based) alternate allele
   **/
  std::string get_alt(const int i) const { return alts[i - 1]; }

  /**
   * @return the variant quality
   **/
  std::string get_quality() const { return quality; }

  /**
   * Extract the information value associated to a key.
   *
   * @param k is the key
   * @return the value associated to key k
   **/
  std::string get_info(const std::string &k) const { return info.at(k); }
};

/*---------- VCF class ----------*/
//! Class for vcf file
/*!
  VCF file.
*/
class VCF {
private:
  /// htslib file
  htsFile *vcf;
  /// htslib header
  bcf_hdr_t *vcf_header;
  /// htslib record
  bcf1_t *vcf_record;

  /// Total number of samples per variant
  int n_samples;
  /// Variants read
  std::vector<Variant> variants;
  /// Genotype fields as cstring (from "GT" to the end of each line)
  std::vector<str_w_l> fmt_lines;
  /// Number of variants to read at each iteration
  size_t block_size;
  /// Number of variants read in the last iteration
  size_t to_parse;

  /// Number of threads to use
  int n_threads;

public:
  /**
   * Constructor
   *
   * @param vcf_path is the path to the input vcf file
   * @param n_threads_ is the number of threads to use
   * @param block_size_ is the number of variants to read at each iteration
   **/
  VCF(char *vcf_path, const int n_threads_, const int block_size_) {
    vcf = bcf_open(vcf_path, "r");
    vcf_header = bcf_hdr_read(vcf);
    vcf_record = bcf_init();
    vcf_record->max_unpack = BCF_UN_INFO;

    n_samples = bcf_hdr_nsamples(vcf_header);

    n_threads = n_threads_;

    block_size = block_size_;
    to_parse = 0;

    // Allocate space for fmt_lines
    for (size_t i = 0; i < block_size; ++i) {
      str_w_l s;
      s.size = 1024 * 10;
      s.seq = (char *)calloc(s.size, sizeof(char));

      fmt_lines.push_back(s);
    }
  }
  /**
   * Destructor
   **/
  ~VCF() {
    for (size_t i = 0; i < block_size; ++i) {
      free(fmt_lines[i].seq);
    }
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf);
  }

  /**
   * Parse block_size lines from the input vcf file
   *
   * @return false if less than block_size lines have been read (due
   * to EOF), true otherwise.
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

      // 2. we store the format line
      char *samples = get_samples(vcf->line.s, vcf->line.l + 1);
      size_t samples_size = vcf->line.l - (samples - vcf->line.s);
      while (samples_size > fmt_lines[i].size) {
	cerr << "Reallocating memory for variant " << variants.back().get_idx() << endl;
        fmt_lines[i].size *= 1.6;
        fmt_lines[i].seq = (char *)realloc(fmt_lines[i].seq, fmt_lines[i].size);
      }
      strncpy(fmt_lines[i].seq, samples, samples_size);
      fmt_lines[i].seq[samples_size] = '\0';
      ++i;
    }
    to_parse = i;

    // 3. we parse and store the genotypes in parallel
    fill_genotypes();

    return i == n;
  }

  /**
   * @return the variants read in the last iteration
   **/
  std::vector<Variant> get_variants() const { return variants; }

private:
  /**
   * Extract the genotype information from a vcf line
   *
   * @param fmt_line is a line of the input vcf file
   * @param len is the length of fmt_line
   * @return a cstring representing the genotypes information
   **/
  char *get_samples(char *const fmt_line, const size_t len) {
    /*
     * We need len since fmt_line contains \0 to separate the first 8
     * fields (htslib replaces \t with \0)
     */
    char *gt_start = fmt_line;
    char *const gt_end = fmt_line + len;

    /*
     * From VCF format specification: "The first sub-field must always
     * be the genotype (GT) if it is present. There are no required
     * sub-fields."
     */
    // FIXME: manage :
    for (; gt_start != gt_end && !(*gt_start == '\0' && *(gt_start+1) == 'G' && *(gt_start+2) == 'T' && *(gt_start+3) == '\t');
	 ++gt_start)
      ;
    ++gt_start; // strip the first \0 left by the for

    // FIXME: manage no GT (ie gt_start == gt_end)
    return gt_start;
  }

  /**
   * Extracts information from a cstring representing a genotype
   *
   * @param start is a cstring
   * @return a tuple containing the first allele, the second allele, and the
   *phased flag
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

  /**
   * Extracts and stores the genotype fields associated to the
   * {var_idx}-th variant read in the last iteration
   *
   * @param var_idx is an index
   **/
  void fill_variant(const uint32_t var_idx) {
    assert(var_idx < to_parse);
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

  /**
   * Extracts and store the genotypes of all the variants read in the
   * last iteration
   **/
  void fill_genotypes() {
#pragma omp parallel for num_threads(n_threads) shared(fmt_lines, variants)
    for (uint i = 0; i < to_parse; ++i) {
      fill_variant(i);
    }
  }
};

#endif
