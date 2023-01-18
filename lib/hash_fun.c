#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// BKDR Hash Function
unsigned int BKDRHash(char *str)
{
    unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
    unsigned int hash = 0;
 
    while (*str)
    {
        hash = hash * seed + (*str++);
    }
 
    return (hash & 0x7FFFFFFF);
}
 


static uint16_t
hash_func(const char *string, const uint32_t max_hash_num)
{
    uint16_t hash_res = 0xffff;
    uint16_t hash_num = 0xffff;
    if (max_hash_num > 0xffff) {
        fprintf(stderr, "max hash num great than MAX_uint16\n");
        return 0;
    }
    uint32_t str_len = strlen(string), i = 0;
    uint32_t limit = str_len - 1; 
    uint32_t j = 0;
    for (i = 0; i < limit; i++) {
        hash_num = 0xffff;
        if (string[i]) {
            hash_num ^= string[i];
            hash_num <<= 8;
            hash_num ^= string[++i];
        }
        hash_num <<= j % 2;    
        printf("%hu\n", hash_num);
        if (j % 3 == 0) {
            hash_res ^= hash_num;
        } else if (j % 3 == 1) {
            hash_res |= hash_num;
        } else {
            hash_res ^= hash_num;
        }
        printf("%hu\n", hash_res);
        j++;
    }
    for (;i < str_len; i++) {
        hash_num = 0xffff;
        if (string[i]) {
            hash_num ^= string[i];
            hash_num <<= 8;
            hash_num ^= string[++i];
        }
        hash_res ^= hash_num;
    }
    
    hash_res %= max_hash_num;
    return hash_res;
}


int
main(void)
{
    char s[] = "HHHH";
    uint16_t hn = hash_fun(s, 360);
    printf("%hu\n", hn);
}
