#ifndef Mersenne_Twister_RAND
#define Mersenne_Twister_RAND

/* This is a wrapper over mt19937ar.c to make it act as a library */

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);   /* [0,0xffffffff]-interval */
long genrand_int31(void);   /* [0,0x7fffffff]-interval */
double genrand_real1(void); /* [0,1]-real-interval */
double genrand_real2(void); /* [0,1)-real-interval */
double genrand_real3(void); /* (0,1)-real-interval */
double genrand_res53(void); /* [0,1) with 53-bit resolution */

#define mt_rand genrand_int31
#define mt_srand init_genrand

#endif
