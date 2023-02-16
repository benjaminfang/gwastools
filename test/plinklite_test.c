#include "../lib/plinklite.h"

int
main(int argc, char *argv[])
{
    const char *plink_fname = argv[1];
    PLINKFILE plink_data = plinkopen(plink_fname);
    unsigned int indi_num = plink_data.individual_num;
    unsigned int vari_num = plink_data.variant_num;
/*    
    printf("indi_num: %u\n", indi_num);
    printf("vari_num: %u\n", vari_num);

    FAM_LINE fam_line;
    famreadline(plink_data, &fam_line);
    printf("%s %s\n",fam_line.family_id, fam_line.within_famid);

    FAM_LINE_ptr all_fam_line = (FAM_LINE_ptr)malloc(sizeof(FAM_LINE) * indi_num);
    famreadlines(plink_data, all_fam_line, indi_num);
    for (int i = 0; i < indi_num; i++) {
        printf("%s %s %s %s %d %d\n", all_fam_line[i].family_id,
            all_fam_line[i].within_famid,
            all_fam_line[i].father_id, all_fam_line[i].mother_id,
            all_fam_line[i].sex, all_fam_line[i].phenotype_value);
        
    }

    BIM_LINE bim_line;
    bimreadline(plink_data, &bim_line);
    printf("%d %s %s\n", bim_line.chrom, bim_line.rsid, bim_line.allel1);
    BIM_LINE_ptr all_bim_line = (BIM_LINE_ptr)malloc(sizeof(BIM_LINE) * vari_num);
    bimreadlines(plink_data, all_bim_line, vari_num);

    for (int i = 0; i < vari_num; i++) {
        printf("%d %s %s\n", all_bim_line[i].chrom, all_bim_line[i].rsid, all_bim_line[i].allel1);
    }
*/


/*
    char *bed_data = (char *)malloc(sizeof(char) * indi_num);
    for (int j = 0; j < 100; j++) {
        bedreaddata(plink_data, bed_data, indi_num);
        for (int i = 0; i < indi_num; i++) {
            printf("%d ", bed_data[i]);
        }
        printf("\n");
    }

*/

    char *bed_data_5 = (char *)malloc(sizeof(char) * indi_num * 5);
    bedloaddata_n(plink_data, bed_data_5, indi_num * 5, 5, 9);
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < indi_num; j++) {
            printf("%d ", bed_data_5[j]);
        }
        printf("\n");
        bed_data_5 += indi_num;
    }

//    free(all_fam_line);
//    free(all_bim_line);
//    free(bed_data);
    plinkclose(plink_data);
    return 0;
}



