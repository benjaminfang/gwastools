#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../lib/read_epi.h"


int
sort_besd(int argc, char * argv[])
{

    char myname[] = "sort";
    char subcmd[32];
    if (argc > 1 && strcmp() != 0) {
        return 0;
    }


    //parse arguments
    char * filename_prefix = argv[1];

    char epi_filename[512];
    char esi_filename[512];
    char besd_filename[512];

    strncpy(epi_filename, filename_prefix, 500);
    strncpy(esi_filename, filename_prefix, 500);
    strncpy(besd_filename, filename_prefix, 500);
    strcat(epi_filename, ".epi");
    strcat(esi_filename, ".esi");
    strcat(besd_filename, ".besd");

    epi_dt_list * epi_dt = NULL;
    read_epi(epi_filename, &epi_dt);
    
    epi_dt_list ** epi_dt_sorted = NULL;
    int epi_len = sort_epi(epi_dt, &epi_dt_sorted);
 

    free_epi_dt(&epi_dt, &epi_dt_sorted);


    return 1;
}


static int
print_sorted_epi(epi_dt_list ** epi_dt_sorted, int epi_len, const char * fout_filename)
{
    int i = 0;
    char chrom[32];
    char pos[32];
    char ori[8];
    epi_dt_list * dt_node;

    FILE * fout = fopen(fout_filename, "w");

    for (i = 0; i < epi_len; i++) {
        dt_node = epi_dt_sorted[i];
        if (dt_node -> chrom == 0) {
            strcpy(chrom, "NA");
        } else if (dt_node -> chrom == 23) {
            strcpy(chrom, "X");
        } else if (dt_node -> chrom == 24) {
            strcpy(chrom, "Y")
        } else {
            sprintf(chrom, "%u", dt_node -> chrom)
        }

        if (dt_node -> epi_pos == 0) {
            strcpy(pos, "NA");
        } else {
            sprintf(pos, "%u", dt_node -> epi_pos);
        }

        if (dt_node -> ori == 0) {
            strcpy(ori, "NA");
        } else if (dt_node -> ori == 1) {
            strcpy(ori, "+");
        } else if (dt_node -> ori == 2) {
            strcpy(ori, "-");
        }

        fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\n", chrom, dt_node -> epi_id, dt_node -> f3,
            pos, dt_node -> gene_id, ori);

    }

    fclose(fout);
    return 0;
}