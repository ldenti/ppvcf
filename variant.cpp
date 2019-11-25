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

string Variant::get_info(const string &key) const {
  return info_values[info_keys.at(key)];
}

void Variant::store_filter(bcf_hdr_t *header, bcf1_t *record) {
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

void Variant::store_info(bcf_hdr_t *header, bcf1_t *record) {
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
    } else if (type == BCF_BT_FLOAT) {
      float *v = NULL;
      bcf_get_info_float(header, record, key, &v, &ndst);
      value = to_string(*v);
    } else if (type == BCF_BT_CHAR) {
      char *v = NULL;
      bcf_get_info_string(header, record, key, &v, &ndst);
      value = string(v);
    } else {
      cerr << "Unknown type " << type << " (field " << key << ")" << endl;
      exit(1);
    }

    info_keys[key] = info_values.size();
    info_values.push_back(value);

    /** // Maybe we can adapt this
       #define BRANCH(type_t, bcf_ht_t) {				\
       type_t *value = NULL;						\
       bcf_get_info_values(header, record, key, (void**)(&value), &ndst, bcf_ht_t); \
       cout << *value << endl;						\
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
  store_filter(header, record);
  store_info(header, record);
}

void Variant::add_genotype(const GT &gt) {
  assert(gti < nsamples);
  genotypes[gti++] = gt;
}
