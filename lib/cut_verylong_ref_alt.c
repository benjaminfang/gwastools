#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FIELD_BUFFER_LEN 128
#define FIELD_BUFFER_LAST_CHAR_IDX 127
#define FIELD_NUM 7

#include "cut_verylong_ref_alt.h"

typedef struct STRUC {
    int field_num;
    char buffer[FIELD_BUFFER_LEN];
    struct STRUC * next;

} long_line_buffer;


static int print_line_buffer(long_line_buffer *line_buffer, FILE *esi_fout,
                             FILE *field_fout, unsigned long line_num);


int
cut_esi_verylong_field(char *esi_filename, char *esi_filename_cut, char *esi_field_index_filename)
{
    int exit_status = 0;

    FILE * esi_fin = fopen(esi_filename, "r");
    if (! esi_fin) {
        fprintf(stderr, "Error, open %s failed.\n", esi_filename);
        exit_status = 1;
    }

    FILE * esi_fout = fopen(esi_filename_cut, "w");
    if (! esi_fout) {
        fprintf(stderr, "Error, open %s failed.\n", esi_filename_cut);
        exit_status = 1;
    }

    FILE * field_fout = fopen(esi_field_index_filename, "w");
    if (! field_fout) {
        fprintf(stderr, "Error, open %s failed.\n", esi_field_index_filename);
        exit_status = 1;
    }

    if (exit_status) {
        return exit_status;
    }

    long_line_buffer * line_head = (long_line_buffer *)malloc(sizeof(long_line_buffer));
    line_head -> field_num = 0;
    line_head -> next = NULL;
    line_head -> buffer[FIELD_BUFFER_LAST_CHAR_IDX] = '\0';

    long_line_buffer * new_node = NULL;
    long_line_buffer * ptr = NULL;

    ptr = line_head;

    unsigned long line_len = 0;
    unsigned long line_num = 0;
    int field_count = 0;
    int cc = 0;
    int field_len = 0;
    
    while ((cc = fgetc(esi_fin)) != EOF) {
        line_len++;
        if (cc == '\n') {
            line_len = 0;
            field_count = 0;
            field_len = 0;
            line_num ++;
            print_line_buffer(line_head, esi_fout, field_fout, line_num);
            ptr = line_head -> next;
            line_head -> next = NULL;
            while (ptr) {
                new_node = ptr -> next;
                free(ptr);
                ptr = new_node;
            }
            ptr = line_head;
            continue;
        }
        
        if (cc == '\t') {
            ptr -> buffer[field_len] = '\0';
            field_count ++;
            field_len = 0;
            new_node = (long_line_buffer *)malloc(sizeof(long_line_buffer));
            new_node -> field_num = field_count;
            new_node -> next = NULL;
            new_node -> buffer[FIELD_BUFFER_LAST_CHAR_IDX] = '\0';
            ptr -> next = new_node;
            ptr = new_node;
            continue;
            
        }
        
        ptr -> buffer[field_len] = (char)cc; 
        
        field_len++;

        if (field_len == FIELD_BUFFER_LEN) {
            new_node = (long_line_buffer *)malloc(sizeof(long_line_buffer));
            new_node -> field_num = field_count;
            new_node -> next = NULL;
            new_node -> buffer[FIELD_BUFFER_LAST_CHAR_IDX] = '\0';
            field_len = 0;
            ptr -> next = new_node;
            ptr = new_node;
        }
    }

    fclose(esi_fin);
    fclose(esi_fout);
    fclose(field_fout);

    return exit_status;
}


static int
print_line_buffer(long_line_buffer * line_buffer, FILE * esi_fout,
    FILE * field_fout, unsigned long line_num)
{
    int exit_status = 0;
    int str_len = 0;
    char print_buffer[129];
    char esi_lines[FIELD_NUM][128];
    int field_num = 0;
    int i = 0;

    while ( line_buffer ) {
        if (line_buffer -> buffer[FIELD_BUFFER_LAST_CHAR_IDX] == '\0') {
            strcpy(esi_lines[line_buffer -> field_num], line_buffer -> buffer);
            line_buffer = line_buffer -> next;
        } else {
            //the "$FR%u_%u" may could overflow.
            sprintf(print_buffer, "$FR_%u_%u", line_num, line_buffer -> field_num);
            strcpy(esi_lines[line_buffer -> field_num], print_buffer);
            field_num = line_buffer -> field_num;
            fprintf(field_fout, "%s\t", print_buffer);
            while (line_buffer && field_num == line_buffer -> field_num) {
                if (line_buffer -> buffer[FIELD_BUFFER_LAST_CHAR_IDX] == '\0') {
                    strncpy(print_buffer, line_buffer -> buffer, FIELD_BUFFER_LEN);
                    fprintf(field_fout, "%s\n", print_buffer);
                } else {
                    fprintf(field_fout, "%s", line_buffer -> buffer);
                }
                line_buffer = line_buffer -> next;
            }   
        }
    }
    
    for (i = 0; i < FIELD_NUM - 1; i++) {
        fprintf(esi_fout, "%s\t", esi_lines[i]);
    }
    fprintf(esi_fout, "%s\n", esi_lines[FIELD_NUM - 1]);


    return exit_status;
}


