#ifndef READ_BESD_HEAD
#define READ_BESD_HEAD

#ifndef READ_BESD_SRC
extern int extract_besd_epi_dense(const char* besd_filename, uint32_t epi_index,
                           float* beta_array, float* se_array);

extern int read_sparse_offset_data(const char* besd_filename, uint64_t* beta_se_offset);

extern int extract_besd_epi_sparse(const char* besd_filename, uint32_t epi_index,
                            int* esi_index_array, float* beta_array,
                            float* se_array, uint64_t* beta_se_offset);
#endif
#endif