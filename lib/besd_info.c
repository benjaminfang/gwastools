#include "besdtype.h"
#include <stdio.h>


int
get_besd_info(char * besd_filename, besd_info * besd_info_dt)
{
    int exit_status = 0;

    besd_info_dt -> file_type = -9;
    besd_info_dt -> sample_size = -9;
    besd_info_dt -> esi_num = -9;
    besd_info_dt -> epi_num = -9;

    int buffer[4];

    FILE * besd_fin = fopen(besd_filename, "r")
    if (!besd_fin) {
        fprintf(stderr, "Error, open %s failed.\n", besd_filename);
        exit_status = 1;
    }

    if (fread(buffer, sizeof(int), 4, besd_fin) != 4) {
        fprintf(stderr, "Error, read first 4 int failed.\n");
        exit_status = 1;
    }
    fclose(besd_fin);

    if (! exit_status) {
        besd_info_dt -> file_type = buffer[0];
        besd_info_dt -> sample_size = buffer[1];
        besd_info_dt -> esi_num = buffer[2];
        besd_info_dt -> epi_num = buffer[3];
    }
    
    return exit_status;
}