#define BESDFILE_SRC

#include "besdfile.h"

int
besdfileopen(const char *fname, BESDFILE_ptr besd_data)
{
    char fname_full[1024];
    int fname_len = strlen(fname);

    FILE *fin = NULL;
    strcpy(fname_full, fname);
    strcat(fname_full, ".epi");
    fin = fopen(fname_full, "r");
    if (!fin) {
        fprintf(stderr, "open %s failed.\n", fname_full);
        return BESD_OPEN_FILE_FAIL;
    }
    besd_data->epi_file = fin;

    fname_full[fname_len] = '\0';
    strcat(fname_full, ".esi");
    fin = fopen(fname_full, "r");
    if (!fin) {
        fprintf(stderr, "open %s failed\n", fname_full);
        return BESD_OPEN_FILE_FAIL;
    }
    besd_data->esi_file = fin;

    fname_full[fname_len] = '\0';
    strcat(fname_full, ".besd");
    fin = fopen(fname_full, "r");
    if (!fin) {
        fprintf(stderr, "open %s failed\n", fname_full);
    }
    besd_data->besd_file = fin;

    // get  epi line nuber.
    fin = besd_data->epi_file;
    uint32_t line_counter = 0;
    int last_char = -1;
    int cc = 0;
    while ((cc = fgetc(fin)) != EOF) {
        if (cc == '\n') {
            line_counter++;
        }
        last_char = cc;
    }
    if (last_char != '\n') {
        fprintf(stderr, "Warning, file not end by new line character.\n");
        line_counter++;
    }
    if (line_counter == 0) {
        fprintf(stderr, "epi file is empty.\n");
        return BESD_FILE_EMPTY;
    }
    besd_data->probe_num = line_counter;
    rewind(fin);

    // get  esi line nuber.
    fin = besd_data->esi_file;
    line_counter = 0;
    last_char = -1;
    cc = 0;
    while ((cc = fgetc(fin)) != EOF) {
        if (cc == '\n') {
            line_counter++;
        }
        last_char = cc;
    }
    if (last_char != '\n') {
        fprintf(stderr, "Warning, file not end by new line character.\n");
        line_counter++;
    }
    if (line_counter == 0) {
        fprintf(stderr, "esi file is empty.\n");
        return BESD_FILE_EMPTY;
    }
    besd_data->variant_num = line_counter;
    rewind(fin);

    //get besd meta
    fin = besd_data->besd_file;
    int buf[16];
    if (fread(buf, sizeof(int), 16, fin) != 16) {
        fprintf(stderr, "read besd file failed.\n");
        return BESD_FILE_READ_FAIL;
    }

    besd_data->file_type = buf[0];
    //printf("file type: %d\n", besd_data->file_type);
    besd_data->sample_num = buf[1];
    int esi_num = buf[2];
    int epi_num = buf[3];
    if (esi_num != besd_data->variant_num) {
        fprintf(stderr, "esi number not consistent. %u %u\n", esi_num,
            besd_data->variant_num);
        return BESD_NUMBER_NOT_MATCH;
    }
    if (epi_num != besd_data->probe_num) {
        fprintf(stderr, "epi number not consistent.\n");
        return BESD_NUMBER_NOT_MATCH;
    }
    if (besd_data->file_type == BESD_FILE_TYPE_SPARSE) {
        uint64_t read_buf[2];
        if (fread(read_buf, sizeof(uint64_t), 2, fin) != 2) {
            fprintf(stderr, "read besd failed.\n");
            return BESD_FILE_READ_FAIL;
        }
        besd_data->value_num = read_buf[0];
        
        besd_data->offset_seek_head = ftell(fin);

        for (uint32_t i = 0; i < epi_num; i++) {
            if (fread(read_buf, sizeof(uint64_t), 2, fin) != 2) {
                fprintf(stderr, "read besd file failed.\n");
                return BESD_FILE_READ_FAIL;
            }
        }
        if (read_buf[1] != besd_data->value_num) {
            fprintf(stderr, "value number not equal last offset value.\n");
            return BESD_NUMBER_NOT_MATCH;
        }

        uint64_t file_size =
            20 * 4 + epi_num * 2 * 8 + besd_data->value_num * 4 * 2;

        fseek(fin, 0, SEEK_END);
        uint64_t file_size_ftell = ftell(fin);
        if (file_size != file_size_ftell) {
            fprintf(stderr, "file size error.\n");
            return BESD_NUMBER_NOT_MATCH;
        }
        besd_data->besd_file_size = file_size;

        fseek(fin, 20 * 4, SEEK_SET);
        besd_data->index_seek_head = 20 * 4 + epi_num * 8 * 2;
        besd_data->beta_se_seek_head =
            20 * 4 + epi_num * 8 * 2 + besd_data->value_num * 4;

        besd_data->current_probe_variant_len = 0;

    } else if (besd_data->file_type == BESD_FILE_TYPE_DENSE) {

    } else {
        fprintf(stderr, "file type not recognized.\n");
        return BESD_FAIL;
    }
   
    return BESD_SUCCESS;
}


void
besdfilerewind(void *besd_data)
{

    return;
}


int
besdfileseek(void *besd_data)
{

    return 0;
}


int
epireadline(void)
{
    return 0;
}


int
epireadlines(void)
{
    return 0;
}


int
esireadline(void)
{
    return 0;
}


int
esireadlines(void)
{
    return 0;
}

int
besdreaddata(BESDFILE_ptr besd_data, uint32_t *index_buf, float *beta_buf,
    float *se_buf, uint32_t buf_len, uint32_t *read_len)
{
    FILE *fin = besd_data->besd_file;
    uint64_t offset_buf[2];
    if (fread(offset_buf, sizeof(uint64_t), 2, fin) != 2) {
        fprintf(stderr, "read besd failed.\n");
        return BESD_FILE_READ_FAIL;
    }
    besd_data->offset_seek_head += 8 * 2;
    uint32_t data_len = offset_buf[1] - offset_buf[0];
    besd_data->current_probe_variant_len = data_len;
    if (buf_len < besd_data->current_probe_variant_len) {
        fprintf(stderr, "alloced buf not enough.\n");
        return BESD_FAIL;
    }

    fseek(fin, besd_data->index_seek_head, SEEK_SET);
    if (fread(index_buf, sizeof(uint32_t), data_len, fin) != data_len) {
        fprintf(stderr, "besd read failed.\n");
        return BESD_FILE_READ_FAIL;
    }
    besd_data->index_seek_head += data_len * 4 * 2;

    size_t read_status = 0;
    fseek(fin, besd_data->beta_se_seek_head, SEEK_SET);
    if (fread(beta_buf, sizeof(float), data_len, fin) != data_len) {
        fprintf(stderr, "besd read failed.\n");
        return BESD_FILE_READ_FAIL;
    }
    if ((read_status = fread(se_buf, sizeof(float), data_len, fin)) != data_len) {
        if (read_status == 0 &&
            (besd_data->beta_se_seek_head == besd_data->besd_file_size)) {
            return BESD_FILE_EOF;
        } else {
            fprintf(stderr, "besd read failed.\n");
            return BESD_FILE_READ_FAIL;
        }

    }
    besd_data->beta_se_seek_head += data_len * 4 * 2;
    fseek(fin, besd_data->offset_seek_head, SEEK_SET);
    *read_len = data_len;

    return BESD_SUCCESS;
}



int
besd_sparse_write_meta(int file_format, int32_t sample_size, uint32_t variant_num,
    uint32_t probe_num, uint64_t vaule_num, uint64_t *offset, FILE *fout)
{
    int meta16[16];
    meta16[0] = file_format;
    meta16[1] = sample_size;
    meta16[2] = variant_num;
    meta16[3] = probe_num;

    for (int i = 4; i < 16; i++) {
        meta16[i] = -9;
    }

    fwrite(meta16, sizeof(int32_t), 16, fout);
    fwrite(&vaule_num, sizeof(uint64_t), 1, fout);
    fwrite(offset, sizeof(uint64_t), probe_num * 2 + 1, fout);
    return 0;
}


int
besd_sparse_write_variant_index(uint32_t *index, uint32_t index_len, FILE *fout)
{
    fwrite(index, sizeof(uint32_t), index_len, fout);
    fwrite(index, sizeof(uint32_t), index_len, fout);
    return 0;
}


int
besd_sparse_write_beta_se_data(float *beta, float *se, uint32_t data_len, FILE *fout)
{
    fwrite(beta, sizeof(float), data_len, fout);
    fwrite(se, sizeof(float), data_len, fout);
    return 0;
}


int
besd_sparse_write(void)
{

    /*
        method one,
        creat tmp files, and then concatnate them.
    */

    return 0;

}