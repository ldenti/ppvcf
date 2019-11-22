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

VCF::VCF(char *vcf_path, const int &nths) {
  vcf = bcf_open(vcf_path, "r");
  vcf_header = bcf_hdr_read(vcf);
  vcf_record = bcf_init();
  vcf_record->max_unpack = BCF_UN_INFO;

  n_samples = bcf_hdr_nsamples(vcf_header);

  n_threads = nths;
}

VCF::~VCF() {
  bcf_hdr_destroy(vcf_header);
  bcf_destroy(vcf_record);
  bcf_close(vcf);
}

bool VCF::parse(const uint32_t n) {
  variants.clear();
  uint32_t i = 0;
  while (i < n && bcf_read(vcf, vcf_header, vcf_record) == 0) {
    // 1. we unpack till INFO field
    bcf_unpack(vcf_record, BCF_UN_INFO);
    variants.push_back(Variant(n_samples));
    variants.back().update_till_info(vcf_header, vcf_record);

    // 2. we store the line
    // FIXME: we should store only the FORMAT fields
    char *samples = get_samples(vcf->line.s, vcf->line.l + 1);
    size_t samples_s = strlen(samples) + 1;
    char *fmt_line = new char[samples_s];
    memcpy(fmt_line, samples, samples_s);

    fmt_lines.push_back(fmt_line);
    ++i;
  }

  fill_genotypes();

  for (uint i = 0; i < fmt_lines.size(); ++i)
    free(fmt_lines[i]);
  fmt_lines.clear();

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
  for (; gt_start != gt_end && (*gt_start != 'G' || *(gt_start + 1) != 'T');
       ++gt_start)
    ;
  return gt_start;
}

GT VCF::extract_genotype(char *start) {
  if (*start != '.')
    return GT();

  int all1, all2;
  char *start2 = start;

  for (; *start2 != '|' && *start2 != '/'; ++start2)
    ;
  bool phased = *start2 == '|';
  *start2 = '\0';
  ++start2;

  all1 = atoi(start);
  all2 = atoi(start2);

  return GT(all1, all2, phased);
}

void VCF::fill_variant(const uint32_t var_idx) {
  char *samples = fmt_lines[var_idx];

  // Drop "GT" at the beginning of the line
  for (; *samples != '\t'; ++samples)
    ;
  ++samples;

  char *start = samples;
  char *end = start;
  bool last = false;
  while (!last) {
    for (; *end != '\t' && *end != '\0'; ++end)
      ;
    last = *end == '\0';
    *end = '\0';
    GT gt = extract_genotype(start);
    variants[var_idx].add_genotype(gt);
    start = end + 1;
    end = start;
  }
  // TODO: manage ":"
}

void VCF::fill_genotypes() {
#pragma omp parallel for num_threads(n_threads) shared(fmt_lines, variants)
  for (uint i = 0; i < fmt_lines.size(); ++i) {
    fill_variant(i);
  }
}
