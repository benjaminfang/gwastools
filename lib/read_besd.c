#define READ_BEAD_SRC 1

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "read_besd.h"
#include "besd_info.h"


/*
    Arguments:
        const char * besd_filename: file name of besd file.
        uint32_t epi_idex: index of epi, the first one is 0.

    Retrun:
        int: return 0 if success, else return 1.
*/
int
extract_besd_epi_dense(const char * besd_filename, uint32_t epi_index, float * beta_array,
    float * se_array)
{
    besd_info besd_info_dt;
    get_besd_info(besd_filename, &besd_info_dt);
    if (besd_info_dt.file_type != DENSE) {
        fprintf(stderr, "File format is not DENSE 5\n");
        return 1;
    }

    int epi_num = besd_info_dt.epi_num;
    int esi_num = besd_info_dt.esi_num;

    FILE * fin = fopen(besd_filename, "r");
    
    long seek_len = 0;
    seek_len = sizeof(int) * 16 + sizeof(float) * esi_num * 2 * epi_index;
    fseek(fin, seek_len, SEEK_SET);
    if (fread(beta_array, sizeof(float), esi_num, fin) != esi_num) {
        fprintf(stderr, "Read besd file failed.\n");
    }
    if (fread(se_array, sizeof(float), esi_num, fin) != esi_num) {
        fprintf(stderr, "Read besd file failed.\n");
    }
    fclose(fin);

    return 0;
}


int
read_sparse_offset_data(const char * besd_filename, uint64_t * beta_se_offset)
{
    besd_info besd_info_dt;
    get_besd_info(besd_filename, &besd_info_dt);
    if (besd_info_dt.file_type != SPARSE) {
      fprintf(stderr, "File format is not SPARSE 3.\n");
    }

    int epi_num = besd_info_dt.epi_num;
    long seek_len = 0;
    seek_len = sizeof(int) * 16 + sizeof(uint64_t) * 1;
    FILE * fin = fopen(besd_filename, "r");
    fseek(fin, seek_len, SEEK_SET);
    int read_len = epi_num * 2 + 1;
    if (fread(beta_se_offset, sizeof(uint64_t), read_len, fin) != read_len) {
        return 1;
    }
    fclose(fin);    
    return 0;
}


int
extract_besd_epi_sparse(const char * besd_filename, uint32_t epi_index, int * esi_index_array,
    float * beta_array, float * se_array, uint64_t * beta_se_offset)
{
    besd_info besd_info_dt;
    get_besd_info(besd_filename, &besd_info_dt);
    if (besd_info_dt.file_type != SPARSE) {
        fprintf(stderr, "File format is not SPARSE 3.\n");
    }

    int epi_num = besd_info_dt.epi_num;
    int esi_num = besd_info_dt.esi_num;
    FILE * fin = fopen(besd_filename, "r");

    long seek_len = 0;
    seek_len = sizeof(int) * 16;
    fseek(fin, seek_len, SEEK_SET);
    uint64_t val_num = 0;
    if (fread(&val_num, sizeof(uint64_t), 1, fin) != 1) {
        return 1;
    }

    if (! beta_se_offset) {
        beta_se_offset = (uint64_t *) malloc(sizeof(uint64_t) * epi_num * 2 + 1);
        if (fread(beta_se_offset, sizeof(uint64_t), epi_num * 2 + 1, fin) != epi_num * 2 + 1) {
            return 1;
        }
    } else {
        seek_len = sizeof(uint64_t) * epi_num * 2 + 1;
        fseek(fin, seek_len, SEEK_CUR);
    }
    
    seek_len = sizeof(int) * beta_se_offset[epi_index * 2];
    int read_len = beta_se_offset[epi_index * 2 + 1] - beta_se_offset[epi_index * 2];
    fseek(fin, seek_len, SEEK_CUR);
    if (fread(esi_index_array, size_t(int), read_len, fin) != read_len) {
        return 1;
    }

    //seek_len = sizeof(int) * (val_num - (beta_se_offset[epi_index * 2] + read_len)) + \
    //    sizeof(int) * beta_se_offset[epi_index * 2];
    seek_len = sizeof(int) * (val_num - read_len);
    fseek(fin, seek_len, SEEK_CUR);
    if (fread(beta_array, sizeof(float), read_len, fin) != read_len) {
        return 1;
    }
    if (fread(se_array, sizeof(float), read_len,fin) != read_len) {
        return 1;
    }

    fclose(fin);
    return 0;
}
