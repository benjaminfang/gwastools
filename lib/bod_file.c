#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>


typedef struct oii_data_node {
    char family_id[64];
    char indiv_id[64];
    char parental_id[64];
    char maternal_id[64];
    char sex;
    struct oii_data_node *next;
} OII_DATA_NODE, *OII_DATA_NODE_ptr;

typedef struct {
    uint32_t line_num;
    OII_DATA_NODE_ptr current_node;
    uint32_t current_node_index;
    OII_DATA_NODE_ptr oii_data;
    OII_DATA_NODE_ptr *oii_data_ptr_array;
} OII_DATA, *OII_DATA_ptr;



typedef struct opi_data_node {
    unsigned char chrom; //0 for NA, 200 for X, 201 for Y, 202 for MT 
    char probe_id[64];
    uint64_t position; // 0 for NA
    char gene_id[64];
    char ori; // 0 for NA, 1 for +, 2 for -
    struct opi_data_node * next;
} OPI_DATA_NODE, *OPI_DATA_NODE_ptr;

typedef struct {
    uint32_t line_num;
    OPI_DATA_NODE_ptr current_node;
    uint32_t current_node_index;
    OPI_DATA_NODE_ptr opi_data;
    OPI_DATA_NODE_ptr *opi_data_ptr_array;
} OPI_DATA, *OPI_DATA_ptr;


typedef struct {
    char value_type;
    char data_type;
    uint32_t indiv_num;
    uint32_t probe_num;
    FILE *filehandle;
    uint32_t current_probe;
    uint64_t current_file_pos;
    uint64_t data_len;
    double *data;
} BOD_DATA, *BOD_DATA_ptr;


static int
clean_oii_data(OII_DATA_ptr data_in) {
    if (data_in->oii_data_ptr_array) {
        free(data_in->oii_data_ptr_array);
        data_in->oii_data_ptr_array = NULL;
    }
    if (data_in->current_node) {
        data_in->current_node = NULL;
    }
    OII_DATA_NODE_ptr last_node = NULL, next_node = NULL;
    last_node = data_in->oii_data;
    data_in->oii_data = NULL;
    while(last_node) {
        next_node = last_node -> next;
        free(last_node);
        last_node = next_node;
    }

    free(data_in);
    return 0;
}


static OII_DATA_ptr
read_oii_file(const char * oii_filename)
{
    OII_DATA_ptr data_out = (OII_DATA_ptr) malloc(sizeof(OII_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc failed.\n");
    }
    data_out->line_num = 0;
    data_out->current_node = NULL;
    data_out->current_node_index = 0;
    data_out->oii_data = NULL;
    data_out->oii_data_ptr_array = NULL;
    FILE * fin = fopen(oii_filename, "r");
    if (!fin) {
        printf("%s\n", oii_filename);
        fprintf(stderr, "open oii file failed.\n");
    }
    
    char family_id[64];
    family_id[63] = '\0';
    char indiv_id[64];
    indiv_id[63] = '\0';
    char parental_id[64];
    parental_id[63] = '\0';
    char maternal_id[64];
    maternal_id[63] = '\0';
    char sex[64];
    sex[63] = '\0';
    OII_DATA_NODE_ptr first_node = NULL, new_node = NULL, node_ptr = NULL;
    char line_buf[512];
    line_buf[511] = '\0';
    while(fgets(line_buf, 512, fin)) {
        if (line_buf[511] != '\0') {
            fprintf(stderr, "read oii line buffer overflow.\n");
            return NULL;
        }
        if (strlen(line_buf) > 262) {
            fprintf(stderr, "oii line field overflow.\n");
            return NULL;
        }

        if (sscanf(line_buf, "%s %s %s %s %s", family_id, indiv_id, parental_id,
            maternal_id, sex) != 5) {
            fprintf(stderr, "sscanf failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }

        new_node = (OII_DATA_NODE_ptr) malloc(sizeof(OII_DATA_NODE));
        if (!new_node) {
            fprintf(stderr, "malloc for oii node failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        new_node->next = NULL;

        if (family_id[63] == '\0') {
            strcpy(new_node->family_id, family_id);
        } else {
            fprintf(stderr, "family id field overflow.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        if (indiv_id[63] == '\0') {
            strcpy(new_node->indiv_id, indiv_id);
        } else {
            fprintf(stderr, "oii individual field overflow.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        if (parental_id[63] == '\0') {
            strcpy(new_node->parental_id, parental_id);
        } else {
            fprintf(stderr, "parental id field failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        if (maternal_id[63] == '\0') {
            strcpy(new_node->maternal_id, maternal_id);
        } else {
            fprintf(stderr, "maternal id field failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        if (sex[63] == '\0') {
            if (strcmp(sex, "NA") == 0) {
                new_node->sex = 0;
            } else if (strlen(sex) == 1 && isdigit(*sex)) {
                new_node->sex = (char) atoi(sex);
            } else {
                fprintf(stderr, "sex field not recgnized.\n");
                clean_oii_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "sex field failed.\n");
            clean_oii_data(data_out);
            return NULL;
        }
        
        if (first_node == NULL) {
            first_node = new_node;
            node_ptr = new_node;
            data_out->line_num ++;
            data_out->current_node=first_node;
            data_out->oii_data = first_node;
        } else {
            data_out->line_num ++;
            node_ptr->next = new_node;
            node_ptr = new_node;
        }
    }
    
    return data_out;
}


static int
clean_opi_data(OPI_DATA_ptr data_in)
{
    if (data_in->opi_data_ptr_array) {
        free(data_in->opi_data_ptr_array);
        data_in->opi_data_ptr_array = NULL;
    }
    data_in->current_node = NULL;
    OPI_DATA_NODE_ptr last_node = NULL, next_node = NULL;
    last_node = data_in->opi_data;
    data_in->opi_data = NULL;
    if (last_node) {
        next_node = last_node->next;
        free(last_node);
        last_node = next_node;
    }
    free(data_in);
    return 0;
}


static OPI_DATA_ptr
read_opi_file(const char * opi_filename)
{
    OPI_DATA_ptr data_out = (OPI_DATA_ptr) malloc(sizeof(OPI_DATA));
    if (!data_out) {
        fprintf(stderr, "malloc for opi data failed.\n");
        return NULL;
    }
    data_out->line_num = 0;
    data_out->current_node = NULL;
    data_out->current_node_index = 0;
    data_out->opi_data = NULL;
    data_out->opi_data_ptr_array = NULL;

    FILE * fin = fopen(opi_filename, "r");
    if (!fin) {
        fprintf(stderr, "open opi file failed.\n");
        return NULL;
    }

    char chrom[16];
    chrom[15] = '\0';
    char probe_id[64];
    probe_id[63] = '\0';
    char position[16];
    position[15] = '\0';
    char gene_id[64];
    gene_id[63] = '\0';
    char oritation[16];
    oritation[15] = '\0';
    OPI_DATA_NODE_ptr first_node = NULL, new_node = NULL, node_ptr = NULL;
    char line_buf[512];
    line_buf[511] = '\0';
    char * char_ptr = NULL;
    int str_len = 0;
    while(fgets(line_buf, 512, fin)) {
        if (line_buf[511] != '\0') {
            fprintf(stderr, "read oii line buffer overflow.\n");
            return NULL;
        }
        if (sscanf(line_buf, "%s %s %s %s %s", chrom, probe_id, position, 
            gene_id, oritation) != 5) {
            fprintf(stderr, "opi field failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }

        new_node = (OPI_DATA_NODE_ptr) malloc(sizeof(OPI_DATA_NODE));
        if (!new_node) {
            fprintf(stderr, "malloc for new opi node failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }
        new_node->next = NULL;
        if (chrom[15] == '\0') {
            str_len = strlen(chrom);
            char_ptr = chrom;
            while (isdigit(*char_ptr)) { 
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                new_node->chrom = (unsigned char) atoi(chrom);
            } else if (strcmp(chrom, "X") == 0) {
                new_node->chrom = 200;
            } else if (strcmp(chrom, "Y") == 0) {
                new_node->chrom = 201;
            } else if (strcmp(chrom, "MT") == 0) {
                new_node->chrom = 202;
            } else if (strcmp(chrom, "NA") == 0){
                new_node->chrom = 0;
            } else {
                fprintf(stderr, "chrom not reganized.\n");
                clean_opi_data(data_out);
                return NULL;
            }
        }
        if (probe_id[63] == '\0') {
            strcpy(new_node->probe_id, probe_id);
        } else {
            fprintf(stderr, "probe id field failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }
        if (position[15] == '\0') {
            char_ptr = position;
            str_len = strlen(position);
            while(isdigit(*char_ptr)) {
                str_len--;
                char_ptr++;
            }
            if (!str_len) {
                new_node->position = (uint64_t) atol(position);
            } else if (strcmp(position, "NA") == 0) {
                new_node->position = 0;
            } else {
                fprintf(stderr, "unknow position field.\n");
                clean_opi_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "opi position failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }
        if (gene_id[63] == '\0') {
            strcpy(new_node->gene_id, gene_id);
        } else {
            fprintf(stderr, "gene id field failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }
        if (oritation[15] == '\0') {
            if (strcmp(oritation, "+") == 0) {
                new_node->ori = 1;
            } else if (strcmp(oritation, "-") == 0) {
                new_node->ori = 2;
            } else if (strcmp(oritation, "NA") == 0) {
                new_node->ori = 0;
            } else {
                fprintf(stderr, "oritation not recognized.\n");
                clean_opi_data(data_out);
                return NULL;
            }
        } else {
            fprintf(stderr, "oritation field failed.\n");
            clean_opi_data(data_out);
            return NULL;
        }

        if (first_node == NULL) {
            first_node = new_node;
            node_ptr = new_node;
            data_out->line_num ++;
            data_out->opi_data = first_node;
            data_out->current_node = first_node;

        } else {
            data_out->line_num ++;
            node_ptr->next = new_node;
            node_ptr = new_node;
        }
    }

    return data_out;

}


static BOD_DATA_ptr
init_bod_data_struct(const char * bod_filename)
{

    BOD_DATA_ptr data_out = (BOD_DATA_ptr) malloc(sizeof(BOD_DATA));
    if (!data_out){
        fprintf(stderr, "alloc mem failed.\n");
        return NULL;
    }
    data_out->filehandle = NULL;
    data_out->data = NULL;
    data_out->data_len = 0;

    FILE * fin = fopen(bod_filename, "rb");
    if (!fin){
        fprintf(stderr, "open file failed.\n");
        return NULL;
    }
    char buf[12];
    if (fread(buf, sizeof(char), 12, fin) != 12) {
        fprintf(stderr, "read bod first 12 bytes faild.\n");
        fclose(fin);
        return NULL;
    }

    data_out->value_type = buf[0];
    data_out->data_type = buf[1];
    data_out->indiv_num = ((uint32_t *)buf)[1];
    data_out->probe_num = ((uint32_t *)buf)[2];
    data_out->filehandle = fin;
    data_out->current_probe = 0;
    data_out->current_file_pos = 12;
    data_out->data = (double *)malloc(sizeof(double) * data_out->indiv_num);
    
    return data_out;
}


static int
final_bod_data_struct(BOD_DATA_ptr data_in)
{
    if (data_in->data_len != 0){
        free(data_in->data);
    }
    if (data_in->filehandle) {
        fclose(data_in->filehandle);
    }
    free(data_in);
    return 0;
}


static void
print_bod_info(BOD_DATA_ptr data_in)
{
    printf("value_type:        %d\n", data_in->value_type);
    printf("data_type:         %d\n", data_in->data_type);
    printf("individule number: %u\n", data_in->indiv_num);
    printf("probe number:      %u\n", data_in->probe_num);
    return;
}


static int
read_bod_by_probe(BOD_DATA_ptr data_in)
{
    if (data_in->current_probe >= data_in->probe_num){
        return 1;
    }
    if (fread(data_in->data, sizeof(double), data_in->indiv_num, data_in->filehandle) 
        != data_in->indiv_num){
            fprintf(stderr, "read probe data failed.\n");
            return 2;
    }
    data_in->data_len = data_in->indiv_num;
    data_in->current_probe++;
    data_in->current_file_pos += sizeof(double) * data_in->indiv_num;
    return 0;
}


static BOD_DATA_ptr
read_bod_totall(BOD_DATA_ptr)
{

    return NULL;
}


static BOD_DATA_ptr
query_bod_file(uint32_t probe_id)
{
    return NULL;
}



int
main(int argc, char *argv[])
{
    char file_name[512];
    strcpy(file_name, argv[1]);
    uint32_t file_prefix_len = strlen(file_name);
    strcat(file_name, ".bod");
    BOD_DATA_ptr bod_data = init_bod_data_struct(file_name);
    //print_bod_info(bod_data);
    uint32_t i = 0, j = 0;
    uint32_t ind_num = 0;
    //uint32_t pro_num = 0;
    ind_num = bod_data->indiv_num;
    //pro_num = bod_data->probe_num;
    while (true) {
        int res_status = read_bod_by_probe(bod_data);
    
        for (j = 0; j < ind_num; j++) {
            printf("%lf\n", (bod_data->data)[j]);
        }
        break;
    }
    return 0;
    file_name[file_prefix_len] = '\0';
    strcat(file_name, ".oii");
    OII_DATA_ptr oii_data = read_oii_file(file_name);
    printf("oii line number: %u\n", oii_data->line_num);
    uint32_t oii_line_num = oii_data->line_num;
    OII_DATA_NODE_ptr oii_node_ptr = oii_data->oii_data;
    for (i = 0; i < oii_line_num; i++) {
        continue;
        printf("%s %s %s %s %d\n", oii_node_ptr->family_id, oii_node_ptr->indiv_id,
            oii_node_ptr->parental_id, oii_node_ptr->maternal_id, oii_node_ptr->sex);
        
        oii_node_ptr = oii_node_ptr->next;
    }

    file_name[file_prefix_len] = '\0';
    strcat(file_name, ".opi");
    OPI_DATA_ptr opi_data = read_opi_file(file_name);
    printf("opi line number: %u\n", opi_data->line_num);
    uint32_t opi_line_num = opi_data->line_num;
    OPI_DATA_NODE_ptr opi_node_ptr = opi_data->opi_data;
    for (i = 0; i < opi_line_num; i++) {
        printf("%d %s %lu %s %d\n", opi_node_ptr->chrom, opi_node_ptr->probe_id,
            opi_node_ptr->position, opi_node_ptr->gene_id, opi_node_ptr->ori);
        opi_node_ptr = opi_node_ptr->next;
    }

    final_bod_data_struct(bod_data);
    clean_oii_data(oii_data);
    return 0;
}
