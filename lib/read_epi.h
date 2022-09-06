#ifndef READ_EPI_HEAD
#define READ_EPI_HEAD 1
#include <stdint.h>

typedef struct STRUC
{
    unsigned char chrom; // 0 for NA
    char epi_id[128];
    char f3[32];
    uint32_t epi_pos; // 0 for NA
    char gene_id[128];
    unsigned char ori; // 0 for NA, 1 for +, 2 for -

    struct STRUC *next;

} epi_dt_list;


void * read_epi(const char * epi_filename, epi_dt_list ** epi_dt);
int free_epi_dt(epi_dt_list ** epi_dt);
int sort_epi(epi_dt_list * epi_dt, epi_dt_list *** epi_dt_sorted);

#endif
