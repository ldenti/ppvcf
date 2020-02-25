#ifndef _PPVCF_HPP_
#define _PPVCF_HPP_

#include <assert.h>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string.h>
#include <string>
#include <tuple>
#include <vector>

#include <omp.h>

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
   * Convert the info field in a map<string, string>.
   *
   * @param info is the info field
   **/
  void store_info(const std::string &info_string) {
    std::vector<std::string> info_tokens;
    std::istringstream info_stream(info_string);
    std::string token;
    while (std::getline(info_stream, token, ';'))
      info_tokens.push_back(token);
    for (const auto token : info_tokens) {
      int eq_pos = token.find("=");
      std::string key = token.substr(0, eq_pos);
      std::string value = token.substr(eq_pos + 1);
      info[key] = value;
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
   * Extract and store the first eight fields of a record.
   *
   * @param fields are the 8 mandatory fields of a VCF record
   **/
  void update_till_info(const std::vector<std::string> &fields) {
    seq_name = fields[0];
    pos = stoi(fields[1]) - 1;
    idx = fields[2];
    ref = fields[3];

    std::istringstream alts_stream(fields[4]);
    std::string alt;
    while (std::getline(alts_stream, alt, ',')) {
      alts.push_back(alt);
    }

    quality = fields[5];
    filter = fields[6];

    store_info(fields[7]);
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

  /* Just for testing :) */
  // void print() const {
  //   std::cout << seq_name << "\t"
  // 	      << pos+1 << "\t"
  //     	      << idx << "\t"
  //     	      << ref << "\t";
  //   for(uint i = 0; i<alts.size(); ++i) {
  //     std::cout << alts[i];
  //     if(i != alts.size()-1)
  // 	std::cout << ",";
  //     else
  // 	std::cout << "\t";
  //   }
  //   std::cout << quality << "\t"
  // 	      << filter << "\t";
  //   for (std::map<std::string,std::string>::const_iterator it=info.begin();
  //   it!=info.end(); ++it)
  //     std::cout << it->first << "=" << it->second << ";";
  //   std::cout << "\t" << "GT" << "\t";

  //   for(const auto gt : genotypes)
  //     std::cout << (int)gt.a1 << "|" << (int)gt.a2 << "\t";
  //   std::cout << "\n";
  // }
};

/*---------- VCF class ----------*/
//! Class for vcf file
/*!
  VCF file.
*/
class VCF {
private:
  /// vcf file
  FILE *vcf;
  /// vcf header
  std::string vcf_header;
  /// single vcf record (temp)
  str_w_l line;
  /// vcf records
  std::vector<str_w_l> vcf_records;
  /// Parsed variants
  std::vector<Variant> variants;

  /// Total number of samples per variant
  int n_samples;
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
  VCF(char *vcf_path, const int n_threads_, const int block_size_)
      : block_size(block_size_), to_parse(0), n_threads(n_threads_) {

    vcf = fopen(vcf_path, "r");
    if (vcf == NULL) {
      std::cout << "Error opening VCF file" << std::endl;
      exit(1);
    }

    // Allocate space for single temp line
    line.size = 1024 * 20;
    line.seq = (char *)calloc(line.size, sizeof(char));

    // Allocate space for records
    for (size_t i = 0; i < block_size; ++i) {
      str_w_l s;
      s.size = 1024 * 10;
      s.seq = (char *)calloc(s.size, sizeof(char));
      vcf_records.push_back(s);
    }

    parse_header();
  }
  /**
   * Destructor
   **/
  ~VCF() {
    for (size_t i = 0; i < block_size; ++i) {
      free(vcf_records[i].seq);
    }
    fclose(vcf);
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
    while (i < block_size && fgets(line.seq, line.size, vcf) != NULL) {
      variants.push_back(Variant(n_samples));

      size_t line_len = strlen(line.seq);
      while (line_len > vcf_records[i].size) {
        std::cerr << "Reallocating memory for variant "
                  << variants.back().get_idx() << std::endl;
        vcf_records[i].size *= 1.6;
        vcf_records[i].seq =
            (char *)realloc(vcf_records[i].seq, vcf_records[i].size);
      }
      strncpy(vcf_records[i].seq, line.seq, line_len);
      vcf_records[i].seq[line_len - 1] = '\0'; // -1 to remove trailing \n

      ++i;
    }
    to_parse = i;

    parse_block();

    return i == n;
  }

  /**
   * @return the variants read in the last iteration
   **/
  std::vector<Variant> get_variants() const { return variants; }

private:
  /**
   * Read and store information from VCF header
   **/
  void parse_header() {
    while (fgets(line.seq, line.size, vcf)) {
      vcf_header += line.seq;
      if (line.seq[0] == '#' and line.seq[1] == 'C')
        break;
    }
    if (!(line.seq[0] == '#' and line.seq[1] == 'C')) {
      std::cerr
          << "Malformed VCF file. I cannot find the header line (#CHROM...)"
          << std::endl;
      exit(1);
    }
    int i = 0;
    int n_cols = 0;
    while (line.seq[i] != '\n') {
      n_cols += line.seq[i] == '\t';
      ++i;
    }
    n_samples = n_cols - 8;
  }

  /**
   * Extracts information from a cstring representing a genotype
   *
   * @param start is a cstring
   * @return a tuple containing the first allele, the second allele,
   * and the phased flag
   **/
  std::tuple<uint8_t, uint8_t, bool> extract_genotype(char *start) {
    if (*start == '.')
      return std::tuple<uint8_t, uint8_t, bool>(-1, -1, false);

    uint8_t all1, all2;
    char *start2 = start;

    for (; *start2 != '|' && *start2 != '/'; ++start2)
      ;
    bool phased = *start2 == '|';
    *start2 = '\0';

    all1 = atoi(start);
    all2 = atoi(start2 + 1);

    *start2 = '|'; // don't care, just fix

    return std::move(std::tuple<uint8_t, uint8_t, bool>(all1, all2, phased));
  }

  /**
   * Extracts and stores the genotype fields associated to the
   * {var_idx}-th variant read in the last iteration
   *
   * @param var_idx is an index
   **/
  void fill_variant(const uint32_t var_idx) {
    assert(var_idx < to_parse);

    char *record = vcf_records[var_idx].seq;

    int field = 0;
    char *start = record;
    char *end = start;
    bool last = *end == '\0';

    std::vector<std::string> fixed_fields(8);
    while (field <= 7) {
      for (; *end != '\t' && *end != '\0'; ++end)
        ;

      last = *end == '\0';
      *end = '\0';
      fixed_fields[field] = start;
      ++field;

      *end = last ? '\0' : '\t';
      start = end + 1;
      end = start;
    }

    variants[var_idx].update_till_info(fixed_fields);

    if (last) // No GTs
      return;

    // Drop "GT" at the beginning of the line
    // TODO: manage ":"
    for (; *start != '\t'; ++start)
      ;
    ++start;
    end = start;
    last = *end == '\0';

    std::tuple<uint8_t, uint8_t, bool> gt;
    while (!last) {
      for (; *end != '\t' && *end != '\0'; ++end)
        ;
      last = *end == '\0';
      *end = '\0';
      gt = extract_genotype(start);
      variants[var_idx].add_genotype(std::get<0>(gt), std::get<1>(gt),
                                     std::get<2>(gt));
      *end = last ? '\0' : '\t';
      start = end + 1;
      end = start;
    }
  }

  /**
   * Extracts and store the information associated to each variant of
   * the current block
   **/
  void parse_block() {
#pragma omp parallel for num_threads(n_threads) shared(vcf_records, variants)
    for (uint i = 0; i < to_parse; ++i) {
      fill_variant(i);
    }
  }
};

#endif
