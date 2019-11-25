#include "vcf_file.hpp"

Variant VCF::front() { return variants.front(); }

Variant VCF::back() { return variants.back(); }

vector<Variant>::iterator VCF::begin() { return variants.begin(); }

vector<Variant>::iterator VCF::end() { return variants.end(); }

vector<Variant>::const_iterator VCF::begin() const { return variants.begin(); }

vector<Variant>::const_iterator VCF::end() const { return variants.end(); }

vector<Variant>::const_iterator VCF::cbegin() const {
  return variants.cbegin();
}

vector<Variant>::const_iterator VCF::cend() const { return variants.cend(); }

VCF::VCF(char *vcf_path, const int nths, const int block_size_) {
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
    s.seq = (char *)malloc(s.size);

    fmt_lines.push_back(s);
  }
}

VCF::~VCF() {
  for (size_t i = 0; i < block_size; ++i) {
    free(fmt_lines[i].seq);
  }
  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);
}

bool VCF::parse() {
  variants.clear();
  uint32_t i = 0;
  uint32_t n = block_size;
  while (i < block_size && bcf_read(vcf, vcf_header, vcf_record) == 0) {
    // 1. we unpack till INFO field
    bcf_unpack(vcf_record, BCF_UN_INFO);
    variants.push_back(Variant(n_samples));
    variants.back().update_till_info(vcf_header, vcf_record);

    // 2. we store the line
    // FIXME: we should store only the FORMAT fields
    char *samples = get_samples(vcf->line.s, vcf->line.l + 1);
    size_t samples_s = vcf->line.l - (samples - vcf->line.s);
    while (samples_s > fmt_lines[i].size) {
      fmt_lines[i].size *= 1.6;
      fmt_lines[i].seq = (char *)realloc(fmt_lines[i].seq, fmt_lines[i].size);
    }
    // char *fmt_line = new char[samples_s];
    strncpy(fmt_lines[i].seq, samples, samples_s);

    // fmt_lines.push_back(string(fmt_line));
    ++i;
  }

  to_parse = i;

  fill_genotypes();

  return i == n;
}

char *VCF::get_samples(char *const fmt_line, const size_t len) {
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

tuple<uint8_t, uint8_t, bool> VCF::extract_genotype(char *start) {
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

void VCF::fill_variant(const uint32_t var_idx) {
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

void VCF::fill_genotypes() {
#pragma omp parallel for num_threads(n_threads) shared(fmt_lines, variants)
  for (uint i = 0; i < to_parse; ++i) {
    fill_variant(i);
  }
}
