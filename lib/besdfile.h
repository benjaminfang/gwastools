#ifndef BESDFILE_HEAD
#define BESDFILE_HEAD
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#ifndef BESDFILE_SRC
#define BESDFILE_EXTERN extern
#else
#define BESDFILE_EXTERN
#endif


#define BESD_SUCCESS 0
#define BESD_FAIL 1
#define BESD_OPEN_FILE_FAIL 2
#define BESD_MALLOC_FAIL 3
#define BESD_FILE_EMPTY 4
#define BESD_FILE_READ_FAIL 5
#define BESD_SAMPLE_NUM_NA -9
#define BESD_NUMBER_NOT_MATCH 6
#define BESD_FILE_EOF -1

#define BESD_FILE_TYPE_SPARSE 3
#define BESD_FILE_TYPE_DENSE 5

typedef struct besdfile_stu {
    uint32_t variant_num;
    uint32_t probe_num;

    FILE *epi_file;
    FILE *esi_file;
    FILE *besd_file;

    uint32_t current_epi_line_index;
    uint32_t current_esi_line_index;
    uint32_t current_besd_probe_data_index;
    
    int file_type;
    int sample_num;
    uint64_t value_num;
    
    uint32_t current_probe_variant_len;
    uint64_t offset_seek_head;
    uint64_t index_seek_head;
    uint64_t beta_se_seek_head;
    uint64_t besd_file_size;

} BESDFILE, *BESDFILE_ptr;

BESDFILE_EXTERN int besdfileopen(const char *fname, BESDFILE_ptr besd_data);
BESDFILE_EXTERN int besdreaddata(BESDFILE_ptr besd_data, uint32_t *index_buf,
                                 float *beta_buf, float *se_buf,
                                 uint32_t buf_len, uint32_t *read_len);

BESDFILE_EXTERN int besd_sparse_write_meta(int file_format, int32_t sample_size,
                        uint32_t variant_num,
                        uint32_t probe_num,
                        uint64_t vaule_num, uint64_t *offset,
                        FILE *fout);
BESDFILE_EXTERN int besd_sparse_write_variant_index(uint32_t *index, uint32_t index_len,
                                    FILE *fout);

BESDFILE_EXTERN int besd_sparse_write_beta_se_data(float *beta, float *se,
                                                   uint32_t data_len,
                                                   FILE *fout);

#endif