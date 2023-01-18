

#define CharM1 0x03
#define CharM2 0x0c
#define CharM3 0x30
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

#define CONVERT_GENO(x) ((x) == 0) ? 2 : (((x) == 2) ? 1 : (((x) == 3) ? 0 : 4))

typedef struct bim_node {
    uint8_t chrom;
    char rsid[GENO_ID_MEM_LEN];
    float phy_pos;
    uint32_t pos;
    char allel1[ALLEL_MEM_LEN];
    char allel2[ALLEL_MEM_LEN];
    struct bim_node *next;

} BIM_NODE, *BIM_NODE_ptr;

typedef struct {
    uint32_t line_num;
    uint32_t current_node_index;
    BIM_NODE_ptr current_node;
    BIM_NODE_ptr bim_data;
    BIM_NODE_ptr *bim_data_ptr_array;
} BIM_DATA, *BIM_DATA_ptr;

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
    FAM_NODE_ptr current_node;
    uint32_t current_node_index;
    FAM_NODE_ptr fam_data;
    FAM_NODE_ptr *fam_data_ptr_array;
} FAM_DATA, *FAM_DATA_ptr;

typedef struct {
    char magic_char[3];
    uint32_t indiv_num;
    uint32_t variant_num;
    FILE *filehandle;
    uint32_t current_variant;
    uint64_t current_file_pos;
    uint64_t data_len;
    char *data;
    uint32_t read_buf_len;
    uint8_t *read_buf;
} BED_DATA, *BED_DATA_ptr;
