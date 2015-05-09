#ifndef SIMPLE_D_TAB_SEPARATED_VALUES_READER
#define SIMPLE_D_TAB_SEPARATED_VALUES_READER

#include <stdio.h>

/* remember to also change read_dTSV_file() if ITEM is changed. */
#define ITEM int

ITEM* read_dTSV_file(FILE* fp, int* out_rows, int* out_cols);
ITEM* read_dTSV(char* filename, int* out_rows, int* out_cols);

#endif
