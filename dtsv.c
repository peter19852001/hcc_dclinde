/*******************************************************************
 A simple function to read tab or space separated files with integer numbers,
 where each line has the same number of numbers.

 Copied and modified from tsv.c
*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dtsv.h"

static void* Malloc_buf(size_t s, char* purpose) {
  void* r = NULL;
  r = malloc(s);
  if(r == NULL) {
    fprintf(stderr, "Error: Cannot allocate %d bytes for %s.\n", s, purpose);
    exit(1);
  }
  return r;
}
#define Malloc(n,t,p) ((t *)Malloc_buf(sizeof(t)*(n), p))

#define N_BLOCK_SIZE 1000
typedef struct s_block {
  /* link list of integers */
  struct s_block* next;
  int n;
  ITEM numbers[N_BLOCK_SIZE];
} block;

static block* add_number_to_block(ITEM d, block* lst) {
  /* add a number to the end of the block if there is space, otherwise allocate new block.
     return the head of the (possibly new) linked blocks. */
  block* p=NULL;
  if(lst==NULL || lst->n >= N_BLOCK_SIZE) { /* need new block */
    p = Malloc(1, block, "add_number_to_block()");
    memset(p,0,sizeof(block));
    p->next = lst;
    p->n = 0;
    lst = p;
  }
  lst->numbers[lst->n++] = d;
  return lst;
}
static void free_blocks(block* lst) {
  block* p=NULL;
  while(lst != NULL) {
    p = lst;
    lst = lst->next;
    free(p);
  }
}

static int read_item(FILE* f, char* buf, int buf_len) {
  /* skip space/tab, then read something until space/tab.
     return EOF if end of file reached before anything useful read.
     return '\n' if newline reached.
     return '\0' otherwise.
     The thing read may be longer than buf_len-1, but only the first buf_len-1 will be read,
     and the rest will be discarded. If anything useful is read, buf will be null-terminated.
  */
  int c=0, i=0;
  while(((c=fgetc(f)) != EOF) && (c==' ' || c=='\t' || c=='\r'))
    ;
  if(c==EOF || feof(f)) return EOF;
  if(c=='\n') return '\n';
  buf[i++] = c;
  while((c=fgetc(f))!=EOF && c!=' ' && c!='\t' && c!='\r' && c!='\n') {
    if(i < buf_len-1) buf[i++] = c;
  }
  buf[i] = '\0';
  if(c=='\n') ungetc(c,f);
  return '\0';
}

/******************************************************************/
ITEM* read_dTSV_file(FILE* fp, int* out_rows, int* out_cols) {
  /* The file is some items separated by tab or space, and each line is assumed to have the same number of items.
     Read the items, and output the number of rows, columns in out_rows and out_cols respectively.
     Returns the newly allocated array of items of length rows*cols if OK. Returns NULL if problem.
     The items will be arranged in row major order, i.e. the (i,j)th entry is at j+i*cols.
  */

#define BUF_LEN 64
  char buf[BUF_LEN];
  block* numbers=NULL;
  block* lst=NULL;
  int nline=0, ncols=-1, n_in_line=0, n_records=0, c=0;
  int i=0,j=0;
  ITEM d=0;
  ITEM* res=NULL;

  if(fp == NULL) return NULL;

  while(1) {
    c = read_item(fp, buf, BUF_LEN);
    if(c == '\n' || c == EOF) { /* finished a line */
      if(n_in_line > 0) {
	n_records++;
	if(ncols < 0) { /* first completed line, get the number of terms */
	  ncols = n_in_line;
	} else if(n_in_line != ncols) {
	  fprintf(stderr, "Error: Inconsistent number of terms in line %d, expected %d, got %d.\n",
		  nline+1, ncols, n_in_line);
	  goto error_exit;
	}
      }
      /* allow empty lines */
      nline++;
      n_in_line = 0; /* prepare for next line */
      if(c == EOF) break;
    } else { /* another term on the same line */
      if(sscanf(buf,"%d", &d) != 1) { /* integer for now, remember to change it when changing ITEM */
	fprintf(stderr, "Error: At line %d, expect a number, got %s.", nline+1, buf);
	goto error_exit;
      }
      numbers = add_number_to_block(d, numbers);
      n_in_line++;
    }
  }
  /* up to here, end of file, got the numbers */
  res = Malloc(ncols*n_records, ITEM, "read_dTSV_file()");
  /* go through the numbers in reverse order, fill them in res (backwards) linearly */
  lst = numbers;
  j = ncols*n_records;

  while(lst != NULL && j>0) {
    for(i = lst->n-1; i>=0; i--) {
      res[--j] = lst->numbers[i];
    }
    lst = lst->next;
  }
  /* done */
  free_blocks(numbers);
  *out_rows = n_records;
  *out_cols = ncols;
  return res;

 error_exit:
  free_blocks(numbers);
  return NULL;
}

ITEM* read_dTSV(char* filename, int* out_rows, int* out_cols) {
  /* Wrapper over read_dTSV_file() */
  int* res=NULL;
  FILE *fp = fopen(filename,"r");
  if(fp == NULL) {
    fprintf(stderr, "Error: Cannot open file %s for reading tab separated values.\n", filename);
    return NULL;
  }
  res = read_dTSV_file(fp, out_rows, out_cols);
  fclose(fp);
  return res;
}

/******************************************************************/
