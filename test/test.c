#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../lib/cut_verylong_ref_alt.h"
#define TEST_cut_verylong_ref_alt 1



int test_cut_verylong_ref_alt(int, char *[]);

int
main(int argc, char * argv[])
{

#ifdef TEST_cut_verylong_ref_alt
    test_cut_verylong_ref_alt(argc, argv);
#endif



        return 0;
}

int
test_cut_verylong_ref_alt(int argc, char * argv[])
{
    int test_status = 0;
    char myname[] = "cut_verlong_ref_alt";
    if (strcmp(argv[1], myname) != 0) {
        return test_status;
    }
    if (argc != 5) {
        exit(0);
    }
    printf("processing...\n");
    cut_esi_verylong_field(argv[2], argv[3], argv[4]);

    return test_status;
}
