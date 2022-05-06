#ifndef BESD_INFO_HEAD
#define BESD_INFO_HEAD

#define DENSE 5
#define SMR_DENSE_3 5

#define SPARSE 3
#define SMR_SPARSE_3 3
#define SMR_SPARSE_3F 0x40400000

// next version of besd dense and sparse format
#define DENSE2 13
#define SPARSE2 14

typedef struct {
  int file_type;
  int sample_size;
  int esi_num;
  int epi_num;

} besd_info;

#endif