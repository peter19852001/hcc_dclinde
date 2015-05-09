#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "parse_option.h"

/*********************************
  Created by Peter Lo
  8 June, 2009
*********************************/

/******************************************************************************/
option* find_option(char* op, option* opts) {
  /* return the pointer to the option op in opts, if not found, returns NULL */
  int i;
  for(i=0; opts[i].str != NULL; i++) {
    if(strcmp(op, opts[i].str)==0) return opts+i;
  }
  return NULL;
}

int is_option_set(char* op, option* opts) {
  /* If the option op is in opts and its index is > 0, which means
     it appears in the argument list, returns true, otherwise false. */
  option* p = find_option(op,opts);
  return (p!=NULL)&&(p->index > 0);
}

int get_int_option_value(char* op, option* opts, int def) {
  /* finds the option op in opts. If found, returns the val_int.
     If not found, returns def. No type checking is done.
     This can be used with both INT and FLAG */
  option* p;
  p = find_option(op, opts);
  return p==NULL ? def : p->val_int;
}
char* get_str_option_value(char* op, option* opts, char* def) {
  /* finds the option op in opts. If found, returns the val_str.
     If not found, returns def. No type checking is done. */
  option* p;
  p = find_option(op, opts);
  return p==NULL ? def : p->val_str;
}
double get_real_option_value(char* op, option* opts, double def) {
  /* finds the option op in opts. If found, returns the val_real.
     If not found, returns def. No type checking is done. */
  option* p;
  p = find_option(op, opts);
  return p==NULL ? def : p->val_real;
}
int get_list_option_range(char* op, option* opts, int* s_idx, int* e_idx) {
  /* finds the option op in opts. If found, returns the start and end
     index in s_idx and e_idx (both inclusive) respectively, and
     returns 1. If not found, returns 0.
  */
  option* p;
  p = find_option(op, opts);
  if(p != NULL) {
    if(s_idx) *s_idx = p->index;
    if(e_idx) *e_idx = p->end_index;
    return 1;
  }
  return 0;
}

int parse_options(int argc, char* argv[], option* opts) {
  /* opts is an array containing the options to be matched, with default values
     filled. End of array opts is marked with an entry with NULL str.
     If an option is matched, its value will be updated.
     If multiple occurrence of an option appears, the one appearing later in
     argv will be used. If no error, then returns 0, 1 otherwise. */
  int i,j;
  option* p;
  /* use simple loops, skipping entry 0 of argv[], since it is the program name */
  for(i=1; i<argc; i++) {
    p = find_option(argv[i],opts);
    if(p!=NULL) { /* matched */
      p->index = i;
      switch(p->type) {
      case FLAG: p->val_int = 1; break;
      case LIST: /* read until an item beginning with "-" */
	for(++i; i<argc && argv[i][0]!='-'; i++); /* empty body */
	p->end_index = --i; /* i now points to either an item beginning with "-" or end of argv, decrement it, so that the outer loop can proceed as normal. */
	break;
      case INT:  p->val_int = atoi(argv[++i]); break;
      case REAL: p->val_real = atof(argv[++i]); break;
      default: /* case STRING: */ p->val_str = argv[++i]; break;
      }
    }
  }
  /* now check that all required options are there */
  i = 0;  /* assuming no error first */
  for(j=0; opts[j].str != NULL; j++) {
    if(opts[j].is_required && opts[j].index <= 0) {
      fprintf(stderr, "Error: Option %s required but not given.\n\t%s:\t%s",
	      opts[j].str, opts[j].str, opts[j].description);
      i = 1;  /* error detected. */
    }
  }
  return i;
}

/*****************************************************************************/

static void pr_option_form(char* buf, option* op) {
  /* Prints the string of the option tag, i.e. what would appear on the first
     line of the usage, into buf. Assumes that buf has enough space.
  */
  char* t="";
  switch(op->type) {
  case FLAG: break;
  case LIST: t = " item1 [item2 ...]"; break;
  case INT: t = " int"; break;
  case REAL: t = " real"; break;
  default: /* case STRING: */ t = " str"; break;
  }
  sprintf(buf,"%s%s%s%s",
	  (op->is_required ? "" : "["),
	  op->str,
	  t,
	  (op->is_required ? "" : "]") );
}

#define IDENT_STR "    "
#define IDENT 4
#define COL_MAX 79
#define USAGE_PROMPT "Usage: "

static char* pr_line(FILE* f, char* s, int e) {
  /* Print a line of s to f. The line would not exceed e columns.
     Returns a pointer to the next portion to be printed if there
     are still characters to be printed. Assuming e > 1.
     Returns NULL if there is nothing more to print.
   */
  int i,c, break_large=0;
  /* Try to look for ' ', '\t', '\n' and '\0' */
  for(i=0,c=0; i<e; i++) {
    if(s[i]==' ' || s[i]=='\t' || s[i]=='\n') c=i;
    if(s[i]=='\0') {c=i; break;}
  }
  if(c==0 && s[0]!='\0') {
    /* a very long word, break it and insert hyphen */
    break_large = 1;
    c = e-1;
  }
  /* Now print the first c characters of s */
  for(i=0; i<c; i++) fputc(s[i],f);
  if(break_large) fputc('-',f);
  fputc('\n',f);
  return s[c]=='\0' ? NULL : (break_large ? s+c : s+c+1);
}
static void pr_description(FILE* f, char* op_str, char* s) {
  /* print the description string op_str and s to f. Each line is idented by IDENT columns.
     Each line would normally not exist COL_MAX columns. Space, tabs and newlines
     are treated as delimiters between words for the purpose of breaking the lines.
   */
  int t;
  char* p;
  t = strlen(op_str);
  fprintf(f,"%s",op_str);
  p = pr_line(f,s, COL_MAX-t);
  while(p != NULL) {
    for(t=0; t<IDENT; t++) fputc(' ',f);
    p = pr_line(f,p,COL_MAX-IDENT);
  }
}

void usage(FILE* f, char* prog_name, option* opts) {
  /* Given the prog_name and the options structure, prints the usage.
     Try to make the layout of the usage look nice.
   */
  int i, c=0, t; /* c tracks the current column */
  char buf[256];
  fprintf(f,"\n%s%s", USAGE_PROMPT, prog_name);
  c = strlen(USAGE_PROMPT) + strlen(prog_name);
  /* now prints the options */
  for(i=0; opts[i].str!=NULL; i++) {
    pr_option_form(buf,opts+i);
    t = strlen(buf);
    if((c+1+t) >= COL_MAX) {
      /* Would go out of the column, print on the next line */
      fprintf(f,"\n%s", buf);
      c = t;
    } else {
      fprintf(f," %s", buf);
      c += t+1;
    }
  }
  fprintf(f,"\n\nDescription of the options:\n");
  /* Now the descriptions */
  for(i=0; opts[i].str!=NULL; i++) {
    sprintf(buf,"  %s%s:  ", opts[i].str, opts[i].is_required ? " [REQUIRED]" : "");
    pr_description(f, buf, opts[i].description);
  }
  fprintf(f,"\n\n");
}
/******************************************************************************/
