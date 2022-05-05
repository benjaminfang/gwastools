#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "read_epi.h"


static int compare_epi(epi_dt_list * p1, epi_dt_list * p2);

void *
read_epi(const char * epi_filename, epi_dt_list ** epi_dt)
{
    int i = 0, j = 0, k = 0;

    char chrom[32];
    char epi_id[128];
    char f3[32];
    char epi_pos[32];
    char gene_id[128];
    char ori[32];

    epi_dt_list * root = NULL, * head = NULL, * new_node = NULL; 

    FILE * fin = fopen(epi_filename, "r");

    if (! fin) {
        fprintf(stderr, "open file failed.\n");
        return NULL;
    }

    while (fscanf(fin, "%s %s %s %s %s %s", chrom, epi_id, f3, epi_pos, gene_id, ori) == 6) {
        new_node = (epi_dt_list *) malloc(sizeof(epi_dt_list));
        new_node -> next = NULL;
        
        j = strlen(chrom);
        if (j == 0) {
            new_node -> chrom = 0;
        } else {
            k = 0;
            for (i = 0; i < j; i++) {
                if (! isdigit(chrom[i])) {
                    k = 1;
                    break;
                }
            }
            if (! k) {
                new_node -> chrom = (unsigned char) atol(chrom);
            } else {
                if (strcmp(chrom, "X") == 0) {
                    new_node -> chrom = 23;
                } else if (strcmp(chrom, "Y") == 0) {
                    new_node -> chrom = 24;
                } else {
                    new_node -> chrom = 0;
                }
            }
        }

        strncpy(new_node -> epi_id, epi_id, 127);
        strncpy(new_node -> f3, f3, 31);
        
        j = strlen(epi_pos);
        k = 0;
        if (j == 0) {
            new_node -> epi_pos = 0;
        } else {
            for (i = 0; i < j; i++) {
                if (!isdigit(epi_pos[i])) {
                    k = 1;
                    break;
                }
            }

            if (!k) {
                new_node -> epi_pos = (uint32_t) atol(epi_pos);
            } else {
                new_node -> epi_pos = 0;
            }
        }

        strncpy(new_node -> gene_id, gene_id, 127);

        if (strcmp(ori, "+") == 0) {
            new_node -> ori = 1;
        } else if (strcmp(ori, "-") == 0) {
            new_node -> ori = 2;
        } else {
            new_node -> ori = 0;
        }

        if (root) {
            head -> next = new_node;
            head = new_node;
        } else {
            root = head = new_node;
        }

    }
    fclose(fin);
    
    if (root) {
        *epi_dt = root;
        return root;
    }

    return NULL;
}


int
free_epi_dt(epi_dt_list ** epi_dt)
{
    epi_dt_list * root = NULL, * head = NULL;
    root = *epi_dt;
    *epi_dt = NULL;
    while (root) {
        head = root -> next;
        free(root);
        root = head;
    }   
    return 0;
}


int
sort_epi(epi_dt_list * epi_dt, epi_dt_list *** epi_dt_sorted)
{
    epi_dt_list ** epi_list_adds = NULL;
    epi_dt_list * root = NULL;
    int epi_dt_len = 0;
    root = epi_dt;
    
    while (root) {
        epi_dt_len++;
        root = root -> next;
    }
    // the address array is 1 greater than length of epi. the last address is NULL.
    epi_list_adds = (epi_dt_list **)malloc(sizeof(epi_dt_list *) * epi_dt_len + 1);

    root = epi_dt;
    
    int i = 0;
    while (root) {
        epi_list_adds[i] = root;
        i++;
        root = root -> next;
    }
    epi_list_adds[i] = NULL;

    epi_dt_list * p1;
    epi_dt_list * p2;
    int j = 0;
    for (i = 1; i < epi_dt_len; i++) {
        j = i - 1;
        p1 = epi_list_adds[i];
        p2 = epi_list_adds[j];
        while (j >= 0 && compare_epi(p1, p2) < 0) {
            epi_list_adds[j + 1] = p2;
            p2 = epi_list_adds[--j];
        }
        epi_list_adds[j + 1] = p1;
    }

    *epi_dt_sorted = epi_list_adds;
    
    return epi_dt_len;
}


static int
compare_epi(epi_dt_list * p1, epi_dt_list * p2)
{
    int cmp_res = 0;
    unsigned char p1_chr, p2_chr;
    uint32_t p1_pos, p2_pos;

    p1_chr = p1 -> chrom;
    p2_chr = p2 -> chrom;
    p1_pos = p1 -> epi_pos;
    p2_pos = p2 -> epi_pos;

    if (p1_chr < p2_chr) {
        cmp_res += -2;
    } else if (p1_chr > p2_chr) {
        cmp_res += 2;
    } else {
        cmp_res += 0;
    }

    if (p1_pos < p2_pos) {
        cmp_res += -1;
    } else if (p1_pos > p2_pos) {
        cmp_res += 1;
    } else {
        cmp_res += 0;
    }

    return cmp_res;
}


