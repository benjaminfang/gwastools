#ifndef HASHFUN_HEAD
#define HASHFUN_HEAD

#ifndef HASHFUN_SRC
#define HASHFUN_EXTERN extern
#else
#define HASHFUN_EXTERN 
#endif
HASHFUN_EXTERN unsigned int BKDRHash(char *str);
#endif