#ifndef MY_PARSE_OPTION_UTILITY
#define MY_PARSE_OPTION_UTILITY

/*********************************
  Created by Peter Lo
  8 June, 2009
*********************************/

/***************************************************************/
/* some little facilities for commandline options */
/* the types:
   STRING means taking the next string as value
   FLAG means not take extra value, but act as a boolean, the value is
     stored in val_int, 0 for absent, 1 for present
   LIST means can take one or more items (not beginning with "-"), the list
     ends with an item beginning with "-", or end of the argument list.
   INT means treating the next string as integer, using atoi()
   REAL means treating the next string as double, using atof() */
enum {STRING=0, FLAG, INT, REAL, LIST};
typedef struct _option {
  char* str;  /* the string to match against */
  short type;   /* the type of value wanted, by default is STRING */
  short is_required;  /* if 1, then this option is required; if missing, an error will be raised */
  char* description;  /* to be printed in usage, should contain the trailing '\n' */
  int index;  /* the index into argv[] at which this option appears, should initially be 0,
		 meaning not appeared; if >0, it means this option has appeared. */
  int end_index; /* the end index of the last item (inclusive) of the LIST type. */
  int val_int;    /* should be filled with the default value wanted */
  double val_real;
  char* val_str;  /* to be filled */
} option;

#define get_flag_option_value(x,y,z) get_int_option_value(x,y,z)
option* find_option(char* op, option* opts);
int is_option_set(char* op, option* opts);
int get_int_option_value(char* op, option* opts, int def);
char* get_str_option_value(char* op, option* opts, char* def);
double get_real_option_value(char* op, option* opts, double def);
int get_list_option_range(char* op, option* opts, int* s_idx, int* e_idx);
int parse_options(int argc, char* argv[], option* opts);
void usage(FILE* f, char* prog_name, option* opts);


#endif
