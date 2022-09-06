#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>

#include "read_esi.h"

static int compare_esi(esi_dt_list * p1, esi_dt_list * p2);


void *
read_esi(const char * esi_filename, esi_dt_list ** esi_dt)
{
    int i = 0, j = 0, k = 0;

    char chrom[32];
    char rsid[128];
    char f3[32];
    char rs_pos[64];
    char ref[128];
    char alt[128];
    char f7[32];

    esi_dt_list * root = NULL, * head = NULL, * new_node = NULL;

    FILE * fin = fopen(esi_filename, "r");
    uint32_t node_id = 0;
    while (fscanf(fin, "%s %s %s %s %s %s %s", chrom, rsid, f3, rs_pos, ref, alt, f7) == 7) {
        node_id++;
        new_node = (esi_dt_list *) malloc(sizeof(esi_dt_list));
        new_node -> next = NULL;
        new_node -> id = node_id;

        j = strlen(chrom);
        if (j == 0) {
            new_node -> chrom = 0;
        }
        else {
            k = 0;
            for (i = 0; i < j; i++) {
                if (!isdigit(chrom[i])) {
                    k = 1;
                    break;
                }
            }
            if (!k) {
                new_node -> chrom = (unsigned char)atol(chrom);
            } else {
                if (strcmp(chrom, "X") == 0) {
                    new_node->chrom = 23;
                }
                else if (strcmp(chrom, "Y") == 0) {
                    new_node->chrom = 24;
                }
                else {
                    new_node->chrom = 0;
                }
            }
        }

        strncpy(new_node -> rsid, rsid, 127);
        strncpy(new_node -> f3, f3, 31);

        j = strlen(rs_pos);
        if (j == 0) {
            new_node -> rs_pos = 0;
        } else {
            k = 0;

            for (i = 0; i < j; i++) {
                if (!isdigit(rs_pos[i])) {
                    k = 1;
                    break;
                }
            }
            if (!k) {
                new_node -> rs_pos = (uint32_t)atol(rs_pos);
            } else {
                new_node -> rs_pos = 0;
            }
        }
        strncpy(new_node -> ref, ref, 127);
        strncpy(new_node -> alt, alt, 127);
        strncpy(new_node -> f7, f7, 31);

        if (!root) {
            head -> next = new_node;
            head = new_node;
        } else {
            root = head = new_node;
        }
    }
    fclose(fin);
    
    if (root) {
        *esi_dt = root;
    }
    
    return root;
}



int
free_esi_dt(esi_dt_list ** esi_dt)
{
    esi_dt_list * root = NULL, * head = NULL;

    root = *esi_dt;
    *esi_dt = NULL;
    while (root) {
        head = root -> next;
        free(root);
        root = head;
    }
    return 0;
}


int
sort_esi(esi_dt_list * esi_dt, esi_dt_list *** esi_dt_sorted)
{
    uint32_t esi_dt_len = 0;
    esi_dt_list * root;
    esi_dt_list ** esi_dt_addr;
    uint32_t i = 0, j = 0;

    root = esi_dt;
    while(root) {
        esi_dt_len++;
        root = root -> next;
    }

    esi_dt_addr = (esi_dt_list **) malloc(sizeof(esi_dt_list*) * esi_dt_len + 1);
    root = esi_dt;
 
    while(root) {
        esi_dt_addr[i] = root;
        i++;
        root = root -> next;
    }
    esi_dt_addr[i] = NULL;

    esi_dt_list * p1 = NULL, * p2 = NULL;
    for (i = 1; i < esi_dt_len; i++) {
        p1 = esi_dt_addr[i];
        j = i - 1;
        p2 = esi_dt_addr[j];
        while (j >= 0 && compare_esi(p1, p2) < 0) {
            esi_dt_addr[j + 1] = p2;
            p2 = esi_dt_addr[--j];
        }

        esi_dt_addr[j + 1] = p1;
    }

    *esi_dt_sorted = esi_dt_addr;

    return esi_dt_len;
}


static int
compare_esi(esi_dt_list * p1, esi_dt_list * p2)
{
    int compare_res = 0;

    unsigned char chr_p1 = p1 -> chrom;
    unsigned char chr_p2 = p2 -> chrom;
    
    uint32_t pos_p1 = p1 -> rs_pos;
    uint32_t pos_p2 = p2 -> rs_pos;

    if (chr_p1 < chr_p2) {
        compare_res += -2;
    } else if (chr_p1 > chr_p2 ) {
        compare_res += 2;
    } else {
        compare_res += 0;
    }


    if  (pos_p1 < pos_p2) {
        compare_res += -1;
    } else if (pos_p1 > pos_p2) {
        compare_res += 1;
    } else {
        compare_res += 0;
    }

    return compare_res;
}

