#ifndef READ_BESD_HEAD
#define READ_BESD_HEAD

typedef struct {
  int file_type;
  int sample_size;
  int esi_num;
  int epi_num;

} besd_info;

#ifndef READ_BESD_SRC

int get_besd_info(const char* besd_filename, besd_info* besd_info_dt);

extern int extract_besd_epi_dense(const char* besd_filename, uint32_t epi_index,
                           float* beta_array, float* se_array);

extern int read_sparse_offset_data(const char* besd_filename, uint64_t* beta_se_offset);

extern int extract_besd_epi_sparse(const char* besd_filename, uint32_t epi_index,
                            int* esi_index_array, float* beta_array,
                            float* se_array, uint64_t* beta_se_offset);
#endif
#endif