#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "plink_file.h"

static int
clean_bim_data(BIM_DATA_ptr data_in)
{
    if (!data_in) {
        return 0;
    }
    if (data_in->bim_data_ptr_array) {
        free(data_in->bim_data_ptr_array);
        data_in->bim_data_ptr_array = NULL;
    }
    if (data_in->current_node) {
        data_in->current_node = NULL;
    }
    BIM_NODE_ptr last_node = NULL, next_node = NULL;
    last_node = data_in->bim_data;
    data_in->bim_data = NULL;

    while(last_node) {
        next_node = last_node->next;
        free(last_node);
        last_node = next_node;
        }
    
    free(data_in);
    return 0;
}


static BIM_DATA_ptr
read_bim_file(const char *bim_filename)
{
    BIM_DATA_ptr data_out = (BIM_DATA_ptr) malloc(sizeof(BIM_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc bim data failed.\n");
        return NULL;
    }
    data_out->line_num = 0;
    data_out->current_node_index = 0;
    data_out->current_node = NULL;
    data_out->bim_data = NULL;
    data_out->bim_data_ptr_array = NULL;
    
    FILE * fin = fopen(bim_filename, "r");
    if (!fin) {
        fprintf(stderr, "open bim file failed.\n");
        return NULL;
    }

    char line_buf[2048];
    line_buf[2047] = '\0';
    BIM_NODE_ptr first_node = NULL, new_node = NULL, node_ptr = NULL;
    char chrom[16];
    chrom[15] = '\0';
    char rsid[64];
    rsid[63] = '\0';
    char phy_pos[16];
    phy_pos[15] = '\0';
    char pos[16];
    pos[15] = '\0';
    char allel1[1024];
    allel1[1023] = '\0';
    char allel2[1024];
    allel2[1023] = '\0';
    int str_len = 0;
    char *char_ptr = NULL;
    char *allel1_ptr;
    char *allel2_ptr;
    while (fgets(line_buf, 2048, fin)) {
        if (line_buf[2047] != '\0') {
            fprintf(stderr, "bim line overflow.\n");
            clean_bim_data(data_out);
            return NULL;
        }

        if (sscanf(line_buf, "%s %s %s %s %s %s", chrom, rsid, phy_pos, pos, allel1, allel2) != 6) {
            fprintf(stderr, "split bim field failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        new_node = (BIM_NODE_ptr) malloc(sizeof(BIM_NODE));
        if (!new_node) {
            fprintf(stderr, "malloc for bim node failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        new_node->allel1 = NULL;
        new_node->allel2 = NULL;
        if (chrom[15] == '\0') {
            str_len = strlen(chrom);
            char_ptr = chrom;
            while (isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                new_node->chrom = (unsigned char)atoi(chrom);
            } else if (strcmp(chrom, "X") == 0) {
                new_node->chrom = 200;
            } else if (strcmp(chrom, "Y") == 0) {
                new_node->chrom = 201;
            } else if (strcmp(chrom, "MT") == 0) {
                new_node->chrom = 202;
            } else if (strcmp(chrom, "NA") == 0) {
                new_node->chrom = 0;
            } else {
                fprintf(stderr, "chrom not reganized.\n");
                clean_bim_data(data_out);
                return NULL;
            }
        }

        if (rsid[63] == '\0') {
            strcpy(new_node->rsid, rsid);
        } else {
            fprintf(stderr, "rsid field failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        if (phy_pos[15] == '\0') {
            new_node->phy_pos = (float)atof(phy_pos);
        } else {
            fprintf(stderr, "physical position field failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        if (pos[15] == '\0') {
            str_len = strlen(pos);
            char_ptr = pos;
            while(isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                new_node->pos = (uint32_t)atol(pos);
            } else if (strcmp(pos, "NA") == 0) {
                new_node->pos = 0;
            } else {
                fprintf(stderr, "bim pos not recognized.\n");
                clean_bim_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "bim position field failed.\n");
            clean_bim_data(data_out);
            return NULL;
        }
        if (allel1[1023] == '\0') {
            str_len = strlen(allel1);
            allel1_ptr = (char *)malloc(++str_len);
            if (allel1_ptr) {
                strcpy(allel1_ptr, allel1);
                new_node->allel1 = allel1_ptr;
            } else {
                fprintf(stderr, "malloc for allel1 failed.\n");
                clean_bim_data(data_out);
                return NULL;
            }
        }
        if (allel2[1023] == '\0') {
            str_len = strlen(allel2);
            allel2_ptr = (char *)malloc(++str_len);
            if (allel2_ptr) {
                strcpy(allel2_ptr, allel2);
                new_node->allel2 = allel2_ptr;
            } else {
                fprintf(stderr, "malloc for allel2 failed.\n");
                clean_bim_data(data_out);
                return NULL;
            }
        }

        if (!first_node) {
            first_node = new_node;
            node_ptr = new_node;
            data_out->current_node = new_node;
            data_out->bim_data = new_node;
            data_out->line_num++;
        } else {
            node_ptr->next = new_node;
            node_ptr = new_node;
            data_out->line_num++;
        }
    }
    return data_out;
}


static int
clean_fam_data(FAM_DATA_ptr data_in)
{
    if (!data_in) {
        return 0;
    }
    if (data_in->fam_data_ptr_array) {
        free(data_in->fam_data_ptr_array);
        data_in->fam_data_ptr_array = NULL;
    }
    if (data_in->current_node) {
        data_in->current_node = NULL;
    }
    FAM_NODE_ptr last_node = NULL, next_node = NULL;
    last_node = data_in->fam_data;
    data_in->fam_data = NULL;
    while (last_node) {
        next_node = last_node->next;
        free(last_node);
        last_node = next_node;
    }
    free(data_in);
    return 0;
}


static FAM_DATA_ptr
read_fam_file(const char *fam_filename)
{
    FAM_DATA_ptr data_out = (FAM_DATA_ptr)malloc(sizeof(FAM_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc for fam data failed.\n");
        return NULL;
    }
    data_out->line_num = 0;
    data_out->current_node = NULL;
    data_out->current_node_index = 0;
    data_out->fam_data = NULL;
    data_out->fam_data_ptr_array = NULL;

    FILE *fin = fopen(fam_filename, "r");
    if (!fin) {
        fprintf(stderr, "open fam file failed.\n");
        return NULL;
    }

    char line_buf[512];
    line_buf[511] = '\0';
    char fam_id[64];
    fam_id[63] = '\0';
    char within_famid[64];
    within_famid[63] = '\0';
    char father_id[64];
    father_id[63] = '\0';
    char mother_id[64];
    mother_id[63] = '\0';
    char sex[16];
    sex[15] = '\0';
    char pheno_val[16];
    pheno_val[15] = '\0';
    FAM_NODE_ptr first_node = NULL, new_node = NULL, node_ptr = NULL;

    while (fgets(line_buf, 512, fin)) {
        if (line_buf[511] != '\0') {
            fprintf(stderr, "fam line buffer overflow.\n");
            clean_fam_data(data_out);
            return NULL;
        }

        if (sscanf(line_buf, "%s %s %s %s %s %s") != 6) {
            fprintf(stderr, "split fam line field failed.\n");
            clean_fam_data(data_out);
            return NULL;   
        }
        new_node = (FAM_NODE_ptr)malloc(sizeof(FAM_NODE));
        if (fam_id[63] == '\0') {
            strcpy(new_node->family_id, fam_id);
        } else {
            fprintf(stderr, "family id field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }
        if (within_famid[63] == '\0') {
            strcpy(new_node->within_famid, within_famid);
        } else {
            fprintf(stderr, "within family id field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }
        if (father_id[63] == '\0') {
            strcpy(new_node->father_id, father_id);
        } else {
            fprintf(stderr, "father id field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }
        if (mother_id[63] == '\0') {
            strcpy(new_node->mother_id, mother_id);
        } else {
            fprintf(stderr, "mother id field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }
        if (sex[15] == '\0') {
            if ((strlen(sex) == 1) && (isdigit(sex[0]))) {
                new_node->sex = (char)atoi(sex);
            } else if (strcmp(sex, "NA") == 0) {
                new_node->sex = 0;
            } else {
                fprintf(stderr, "sex field is not recognized.\n");
                clean_fam_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "sex field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }

        if (pheno_val[15] == '\0') {
            if ((strlen(pheno_val) == 1) && (isdigit(pheno_val[0]))) {
                new_node->phenotype_value = (char)atoi(pheno_val);
            } else if (strcmp(pheno_val, "NA") == 0) {
                new_node->phenotype_value = 0;
            } else {
                fprintf(stderr, "phenotype value field failed.\n");
                clean_fam_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "phenotype value field failed.\n");
            clean_fam_data(data_out);
            return NULL;
        }

        if (!first_node) {
            first_node = new_node;
            node_ptr = new_node;
            data_out->current_node = new_node;
            data_out->fam_data = new_node;
            data_out->line_num++;
        } else {
            node_ptr->next = new_node;
            node_ptr = new_node;
            data_out->line_num++;
        }
    }

    return data_out;
}


static BED_DATA_ptr
init_bed_data_struct(const char * bed_filename, uint32_t indiv_num, uint32_t variant_num)
{
    FILE * fin = fopen(bed_filename, "r");
    if (!fin) {
        fprintf(stderr, "open bed file failed.\n");
        return NULL;
    }

    BED_DATA_ptr data_out = (BED_DATA_ptr)malloc(sizeof(BED_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc for BED_DATA failed.\n");
        return NULL;
    }
    data_out->filehandle = NULL;
    data_out->data_len = 0;
    data_out->data = NULL;
    data_out->current_file_pos = 0;
    data_out->current_variant = 0;
    data_out->indiv_num = 0;
    data_out->variant_num = 0;
    data_out->read_buf_len = 0;
    data_out->read_buf = NULL;

    char buf[3];
    if (fread(buf, sizeof(char), 3, fin) != 3) {
        fprintf(stderr, "read bed first 3 char failed.\n");
        fclose(fin);
        return NULL;
    }
    memcpy(data_out->magic_char, buf, 3);
    data_out->indiv_num = indiv_num;
    data_out->variant_num = variant_num;
    data_out->filehandle = fin;
    data_out->current_file_pos = 3;
    data_out->current_variant = 0;
    data_out->data = (char *) malloc(ceil((double)indiv_num / 4) * 4);
    data_out->read_buf_len = ceil((double)indiv_num / 4);
    data_out->read_buf = (uint8_t *) malloc(ceil((double)indiv_num / 4));
    return data_out;
}


static int
read_bed_by_variant(BED_DATA_ptr data_in)
{
    if (data_in->current_variant >= data_in->variant_num) {
        return 1;
    }

    uint32_t read_len = data_in->read_buf_len;
    if (fread(data_in->read_buf, sizeof(uint8_t), data_in->read_buf_len, data_in->filehandle) != data_in->read_buf_len) {
        fprintf(stderr, "read bed file by variant failed.\n");
        return 2;
    }
    data_in->current_file_pos += data_in->read_buf_len;
    data_in->current_variant += 1;

    int i = 0, j = 0;
    uint8_t char_tmp;
    char Low_1_2b = 0,  Low_2_2b = 0, Low_3_2b = 0, Low_4_2b = 0;
    uint8_t *read_buf = data_in->read_buf;
    for (i = 0; i < read_len; i++) {
        char_tmp = read_buf[i];
        Low_1_2b = char_tmp & CharM1;
        Low_2_2b = (char_tmp & CharM2) >> 2;
        Low_3_2b = (char_tmp & CharM3) >> 4;
        Low_4_2b = char_tmp >> 6;
        (data_in->data)[j] = CONVERT_GENO(Low_1_2b);
        (data_in->data)[++j] = CONVERT_GENO(Low_2_2b);
        (data_in->data)[++j] = CONVERT_GENO(Low_3_2b);
        (data_in->data)[++j] = CONVERT_GENO(Low_4_2b);

        j++;
    }
    if (j != data_in->indiv_num) {
        fprintf(stderr, "decode genotype failed.\n");
        return 3;
    }
    data_in->data_len = data_in->indiv_num;
    
    return 0;
}


static int
final_bed_data_struct(BED_DATA_ptr data_in)
{
    if (data_in) {
        if (data_in->filehandle) {
            fclose(data_in->filehandle);
        }
        if (data_in->data) {
            free(data_in->data);
        }
        if (data_in->read_buf) {
            free(data_in->read_buf);
        }
        free(data_in);
    }
    return 0;
}


#define TEST
#ifdef TEST
int
main(int argc, char *argv[])
{
    char file_name[512];
    int file_prefix_len = strlen(argv[1]);
    strcpy(file_name, argv[1]);
    file_name[file_prefix_len] = '\0';
    strcat(file_name, ".bed");

    BED_DATA_ptr data = init_bed_data_struct(file_name, 360, 9596452);
    uint32_t indiv_num = data->indiv_num;
    uint32_t variant_num = data->variant_num;
    int status = read_bed_by_variant(data);

}

#endif
