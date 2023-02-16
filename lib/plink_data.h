#ifndef PLINK_DATA_HEAD
#define PLINK_DATA_HEAD

#ifndef PLINK_DATA_SRC
#define PLINK_EXTERN extern
#endif

#define GENO_ID_MEM_LEN 64
#define GENO_ID_LEN 62
#define MAX_GENO_NUM 4294967295U
#define MAX_PHENO_NUM 4294967295U
#define MAX_INDIVIDUAL_NUM 4294967295U
#define MAX_CHROM_NUM 200
#define CHROM_X 201
#define CHROM_Y 202
#define CHROM_MT 203
#define ALLEL_MEM_LEN 8
#define ALLEL_LEN 6


#define CharM1 0x03
#define CharM2 0x0c
#define CharM3 0x30

#define CONVERT_GENO(x) ((x) == 0) ? 2 : (((x) == 2) ? 1 : (((x) == 3) ? 0 : 4))

#define BIM_DATA_TYPE 1
typedef struct bim_node {
    char chrom;
    char rsid[64];
    float phy_pos;
    uint32_t pos;
    char *allel1;
    char *allel2;
    struct bim_node *next;

} BIM_NODE, *BIM_NODE_ptr;

typedef struct {
    uint32_t line_num;
    uint32_t current_node_index;
    BIM_NODE_ptr next_iter;
    BIM_NODE_ptr bim_data;
    BIM_NODE_ptr *node_ptr_array;
} BIM_DATA, *BIM_DATA_ptr;

#define FAM_DATA_TYPE 2
typedef struct fam_node {
    char family_id[64];
    char within_famid[64];
    char father_id[64];
    char mother_id[64];
    char sex;
    char phenotype_value;
    struct fam_node *next;
} FAM_NODE, *FAM_NODE_ptr;

typedef struct {
    uint32_t line_num;
    FAM_NODE_ptr next_iter;
    uint32_t current_node_index;
    FAM_NODE_ptr fam_data;
    FAM_NODE_ptr *node_ptr_array;
} FAM_DATA, *FAM_DATA_ptr;

typedef struct {
    char magic_char[3];
    uint32_t indiv_num;
    uint32_t variant_num;
    FILE *filehandle;
    uint32_t current_variant;
    uint64_t current_file_pos;
    uint32_t data_len;
    char *data;
    uint32_t read_buf_len;
    uint8_t *read_buf;
} BED_DATA, *BED_DATA_ptr;

BIM_DATA_ptr read_bim_file(const char *bim_filename);
int clean_bim_data(BIM_DATA_ptr data_in);
FAM_DATA_ptr read_fam_file(const char *fam_filename);
int clean_fam_data(FAM_DATA_ptr data_in);
BED_DATA_ptr init_bed_data_struct(const char *bed_filename,
                                         uint32_t indiv_num,
                                         uint32_t variant_num);
char *read_bed_by_variant(BED_DATA_ptr data_in, char *ext_data_array);
int finalize_bed_data_struct(BED_DATA_ptr data_in);

#endif