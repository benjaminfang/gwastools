#include <stdio.h>
#include <stdlib.h>


#include "../lib/besdfile.h"


int
main(int argc, char *argv[])
{
    const char *fname = argv[1];
    BESDFILE besd_data;

    besdfileopen(fname, &besd_data);
    uint32_t variant_num = besd_data.variant_num;
    uint32_t probe_num = besd_data.probe_num;
    uint32_t *index_buf = (uint32_t *)malloc(sizeof(uint32_t) * variant_num);
    float *beta_buf = (float *)malloc(sizeof(float) * variant_num);
    float *se_buf = (float *)malloc(sizeof(float) * variant_num);
    uint32_t read_len;
    besdreaddata(&besd_data, index_buf, beta_buf, se_buf, variant_num, &read_len);
    for (uint32_t i = 0; i < read_len; i++) {
        printf("%u %e %e\n", index_buf[i], beta_buf[i], se_buf[i]);
    }





    return 0;
}
