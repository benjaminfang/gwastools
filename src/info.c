#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../lib/besd.h"
#include "../lib/besd_mis.h"
#include "../lib/read_epi.h"
#include "../lib/read_esi.h"
#include "../lib/read_besd.h"


int
info(int argc, char * argv[])
{
    char myname[] = "info";
    char subcmd[32];
    strcpy(subcmd, argv[1]);

    if (strcmp(subcmd, myname) != 0) {
        return 0;
    }

    char besd_filename_prefix[512];
    char epi_filename[512];
    char esi_filename[512];
    char besd_filename[512];

    strncpy(besd_filename_prefix, argv[2], 500);
    make_besd_filename(besd_filename_prefix, epi_filename, esi_filename, besd_filename);
    




    return 1;
}