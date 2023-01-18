#include <string.h>

int
make_besd_filename(const char * name_prefix, char * epi_filename, char * esi_filename,
    char * besd_filename)
{

    strncpy(epi_filename, name_prefix, 500);
    strcat(epi_filename, ".epi");
    strncpy(esi_filename, name_prefix, 500);
    strcat(esi_filename, ".esi");
    strncpy(besd_filename, name_prefix, 500);
    strcat(besd_filename, ".besd");

    return 0;
}