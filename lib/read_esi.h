typedef struct STRUC {
  uint32_t id;
  unsigned char chrom;  // 0 for NA, X=23, Y=24.
  char rsid[128];
  char f3[32];
  uint32_t rs_pos;  // 0 for NA
  char ref[128];
  char alt[128];
  char f7[32];

  struct STRUC* next;

} esi_dt_list;
