/*
  Use modified CLINDE to learn model similar to dynamic Bayesian network
  time series discrete data, but with hidden common cause(s).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "parse_option.h"
#include "dtsv.h"
#include "mt_rand.h"
/* we include gsl-1.16 in DCLINDE locally for convenience.
   It should be properly configured, built and installed.
 */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>

/* discrete data in form of small integers */
/* the printing functions print_array2d() and print_carray2d() should
   be changed when ELEM is changed */
#define ELEM int

/************************************************************************/
/* mainly for debug */

static void pr_ints(int n, int x[]) {
  int i=0;
  for(i=0; i<n; i++)
    printf(" %d", x[i]);
  printf("\n");
}
static void pr_int_arr(int n, int ns[], int* xs[]) {
  int i=0;
  printf("\n");
  for(i=0; i<n; i++) {
    printf("%d:", i);
    pr_ints(ns[i],xs[i]);
  }
}
static void pr_doubles(int n, double x[]) {
  int i=0;
  for(i=0; i<n; i++)
    printf(" %f", x[i]);
  printf("\n");
}


/***********************************************************************/
static clock_t stop_watch_start = 0;
static void start_stop_watch() {
  printf("Start stop watch.\n");
  stop_watch_start = clock();
}
static void report_stop_watch() {
  printf("Stop watch elapsed time: %f sec\n",
	 (float)(clock() - stop_watch_start)/CLOCKS_PER_SEC);
}

/***********************************************************************/
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

/***********************************************************************/
typedef struct s_array2d {
  int rows,cols;
  ELEM* v;
} array2d;

/* assume row major format */
static ELEM aref(array2d* d, int i, int j) {return d->v[j+i*(d->cols)];}
static void aset(array2d* d, int i, int j, ELEM val) {d->v[j+i*(d->cols)] = val;}

/* assume column major format */
static ELEM caref(array2d* d, int i, int j) {return d->v[i+j*(d->rows)];}
static void caset(array2d* d, int i, int j, ELEM val) {d->v[i+j*(d->rows)] = val;}
static ELEM* get_col(array2d* d, int c) {return d->v+c*(d->rows);}
/****/
static array2d* new_array2d(int n_rows, int n_cols) {
  array2d* r=NULL;
  r = Malloc(1, array2d, "new_array2d()");
  r->rows = n_rows;
  r->cols = n_cols;
  r->v = Malloc(n_rows*n_cols, ELEM, "new_array2d()");
  return r;
}
static void free_array2d(array2d* d) {
  if(!d) return;
  if(d->v) free(d->v);
  free(d);
}

/* these printing functions should be changed when ELEM is changed */
static void print_array2d(FILE* f, array2d* d) {
  int i=0, j=0;
  for(i=0; i<(d->rows); i++) {
    for(j=0; j<(d->cols); j++) {
      /* fprintf(f, "%f\t", aref(d,i,j)); */
      fprintf(f, "%d\t", aref(d,i,j));
    }
    fprintf(f,"\n");
  }
}

static void print_carray2d(FILE* f, array2d* d) {
  int i=0, j=0;
  for(i=0; i<(d->rows); i++) {
    for(j=0; j<(d->cols); j++) {
      /* fprintf(f, "%f\t", caref(d,i,j)); */
      fprintf(f, "%d\t", caref(d,i,j));
    }
    fprintf(f,"\n");
  }
}

/* for passing array2d to GlobalMIT* and GlobalMIT+ */
void extract_carray2d_row(array2d* d, int i, int L, int out[]) {
  /* get row i of d (which is column major) to out, filling at most L elements */
  int j=0;
  for(j=0; j<L && j<(d->cols); j++) {
    out[j] = caref(d,i,j);
  }
}

/***************************************************************/
typedef struct s_link {
  int from, to; /* from --> to */
  int delay;
  double score;
  double test_value;
} links;

/***************************************************************/
/* for GlobalMIT* and GlobalMIT+ to pass the learnt links to use without them knowing struct s_link. */
links* obtain_links(int n_links) {
  if(n_links <= 0) n_links = 1;
  return Malloc(n_links, links, "obtain_links()");
}
void output_one_link(links* Ls, int i, int from, int to, int delay) {
  /* assume Ls is large enough, put the information to slot i. */
  Ls[i].from = from;
  Ls[i].to = to;
  Ls[i].delay = delay;
}

/***************************************************************/
static int min_int_of(int n, int x[]) {
  int i=0;
  int m=x[0];
  for(i=1; i<n; i++)
    if(x[i] < m) m = x[i];
  return m;
}

static int max_int_of(int n, int x[]) {
  int i=0;
  int m=x[0];
  for(i=1; i<n; i++)
    if(x[i] > m) m = x[i];
  return m;
}

static double max_of(int n, double x[]) {
  int i=0;
  double m=x[0];
  for(i=1; i<n; i++)
    if(x[i] > m) m = x[i];
  return m;
}

static int arg_max_of(int n, double x[]) {
  int i=0, r=0;
  double m=x[0];
  for(i=1; i<n; i++)
    if(x[i] > m) {m = x[i]; r = i;}
  return r;
}

/*****************************************************/
/* need temporary buffer for G-test */
static int* tmp_gtest_buf=NULL;
static int L_tmp_gtest_buf=0;

static void free_temp_buffer() {
  if(tmp_gtest_buf) {
    free(tmp_gtest_buf);
    tmp_gtest_buf = NULL;
    L_tmp_gtest_buf = 0;
  }
}

static int* reserve_temp_buffer(int L) {
  /* make sure the buffer is at least of length L.

     The returned pointer is valid up to the next call to
     reserve_temp_buffer() or free_temp_buffer();
   */
  if(L_tmp_gtest_buf >= L) return tmp_gtest_buf;
  free_temp_buffer();
  tmp_gtest_buf = Malloc(L,int,"reserve_temp_buffer()");
  L_tmp_gtest_buf = L;
  return tmp_gtest_buf;
}

/***** G-test,
       refer to http://en.wikipedia.org/wiki/G-test
       work on discrete data, which is in form of small
       non-negative integers
 *****/
void test_func_gtest(int L, int x[], int y[], double* score, double* test_value) {
  double test_stat=0, p_value=0;
  int i=0, u=0,v=0;
  int *cxy=NULL, *cx=NULL, *cy=NULL;
  int mx=0, my=0;
  /* [0, mx): possible states of x. [0, my): possible states of y. */
  mx = 1 + max_int_of(L, x);
  my = 1 + max_int_of(L, y);

  /* the buffers, no need to release here, as it is temporary */
  /* cxy in row major order, i.e. cxy[v + u*my] is count of (x=u, y=v) */
  cx = reserve_temp_buffer(mx+my+mx*my);
  cy = cx + mx;
  cxy = cy + my;
  memset(cx, 0, sizeof(int)*(mx+my+mx*my));

  /* go through x and y to get the counts */
  for(i=0; i<L; i++) {
    u = x[i]; v = y[i];
    cxy[v + u*my]++;
    cx[u]++;
    cy[v]++;
  }
  /* G test statistics, and the real number of states in x and y */
  /* G = 2*\sum_{i}{O_i * ln(O_i/E_i)} for non-empty cells i */
  test_stat = 0;
  i=0;
  for(u=0; u<mx; u++) {
    for(v=0; v<my; v++, i++) {
      if(cxy[i] > 0) {
	/* debug */
	/*
	printf("x: %d\ty: %d\tcount: %d\n", u,v,cxy[i]);
	*/

	test_stat += cxy[i]*log(cxy[i]/(cx[u]*cy[v]/(double)L));
      }
    }
  }
  test_stat *= 2;

  for(u=0, i=0; i<mx; i++) if(cx[i]>0) u++;
  for(v=0, i=0; i<my; i++) if(cy[i]>0) v++;
  /* approximate chi-square distribution, df=(u-1)*(v-1) */
  p_value = gsl_cdf_chisq_Q(test_stat, (u-1)*(v-1));

  /* debug */
  /*
  printf("** test_stat: %g\tp-value: %g\tdf: %d\n", test_stat,p_value,(u-1)*(v-1));
  */

  if(score) *score = -log10(p_value);
  if(test_value) *test_value = test_stat;
}

void test_c_func_gtest(int L, int x[], int y[], int nZ, int Z[],
		      double* score,double* test_value) {
  double test_stat=0, p_value=0;
  int i=0, j=0, u=0,v=0,w=0, r=0;
  int *czxy=NULL, *czx=NULL, *czy=NULL, *cz=NULL;
  int mx=0, my=0, mz=0, *mzs=NULL;
  /* [0, mx): possible states of x. */
  /* [0, my): possible states of y. */
  /* [0, mz): possible states of Z. */
  mx = 1 + max_int_of(L, x);
  my = 1 + max_int_of(L, y);
  mz = 1;
  for(i=0; i<nZ; i++) {
    mz *= (1 + max_int_of(L, Z + i*L));
  }

  /* the buffers, no need to release here, as it is temporary */
  /* czxy in row major order,
     i.e. cxy[v + u*my + w*mx*my] is count of (Z=w, x=u, y=v) */
  /* czx[u + w*mx]  is count of (Z=w, x=u) */
  /* czy[v + w*my]  is count of (Z=w, y=v) */
  i = nZ + mz + mx*mz + my*mz + mx*my*mz;
  mzs = reserve_temp_buffer(i);  /* mzs length nZ */
  memset(mzs, 0, sizeof(int)*i);
  cz = mzs + nZ;                 /* cz length mz */
  czx = cz + mz;                 /* czx length mx*mz*/
  czy = czx + mx*mz;             /* czy length my*mz */
  czxy = czy + my*mz; /* czxy has length mx*my*mz */

  /* get the number of possible states of each z again */
  for(i=0; i<nZ; i++) {
    mzs[i] = 1 + max_int_of(L, Z + i*L);
  }

  /* go through x, y and Z to get the counts */
  for(i=0; i<L; i++) {
    u = x[i]; v = y[i];
    /* get an index for Z */
    w = 0;
    for(r=1, j=0; j<nZ; j++) {
      w += r*Z[i + j*L];
      r *= mzs[j];
    }

    czxy[v + u*my + w*mx*my]++;
    czx[u + w*mx]++;
    czy[v + w*my]++;
    cz[w]++;
  }
  /* G test statistics, and the real number of states in x, y and Z */
  /* G = 2*\sum_{i}{O_i * ln(O_i/E_i)} for non-empty cells i */
  test_stat = 0;

  for(w=0; w<mz; w++) {
    if(cz[w] <= 0) continue; /* skip empty Z states */
    i = w*mx*my;
    for(u=0; u<mx; u++) {
      for(v=0; v<my; v++, i++) {
	if(czxy[i] > 0) {
	  /* debug */
	  /*
	  printf("z: %d\tx: %d\ty: %d\tcount: %d\n", w,u,v,czxy[i]);
	  */

	  test_stat += czxy[i]*log(czxy[i]/(czx[u + w*mx]*czy[v + w*my]/(double)cz[w]));
	}
      }
    }
  }
  test_stat *= 2;

  for(u=0, i=0; i<mx; i++) {
    for(j=0; j<mz; j++)
      if(czx[i + j*mx]>0) {u++; break;}
  }
  for(v=0, i=0; i<my; i++) {
    for(j=0; j<mz; j++)
      if(czy[i + j*my]>0) {v++; break;}
  }
  for(w=0, i=0; i<mz; i++) if(cz[i]>0) w++;
  /* approximate chi-square distribution, df=w*(u-1)*(v-1) */
  p_value = gsl_cdf_chisq_Q(test_stat, w*(u-1)*(v-1));

  /* debug */
  /*
  printf("** test_stat: %g\tp-value: %g\tdf: %d\n", test_stat,p_value,w*(u-1)*(v-1));
  */

  if(score) *score = -log10(p_value);
  if(test_value) *test_value = test_stat;
}

/***************************************************************/
#define M_DCLINDE_NAME  "dclinde"
#define M_MIT_NAME      "mit"
#define M_MIT_FAST_NAME "mit_fast"

enum methods {METHOD_GLOBALMIT_FAST=0, METHOD_GLOBALMIT, METHOD_DCLINDE,  N_METHODS};
char* method_names[N_METHODS] = {M_MIT_FAST_NAME, M_MIT_NAME, M_DCLINDE_NAME};
char* method_long_names[N_METHODS] = {"GlobalMIT*", "GlobalMIT+", "D-CLINDE"};

enum prunings {PRUNING_ALL=0, PRUNING_COMMON,  N_PRUNINGS};
char* pruning_names[N_PRUNINGS] = {"all", "common"};
/***********/

static int shift_by(ELEM x[], int n_segs, int lens[],
		    int lag, int nL, ELEM outx[]) {
  /* each segment of x is shifted by lag, and each segment's
     effective length is subtracted by nL.
     If a segment is too short, it is not included.

     Returns the length of the shifted series.
  */
  int i=0, j=0, k=0, n=0;
  for(j=0, k=0; k<n_segs; j+=lens[k++]) {
    n = lens[k] - nL;
    if(n <= 0) continue; /* segment too short */
    memcpy(outx+i, x+j+lag, sizeof(ELEM)*n);
    i += n;
  }
  return i;
}

static int shift_by_offset(ELEM x[], int n_segs, int lens[],
			   int offset, int eL, int nL, ELEM outx[]) {
  /* each segment of x has extra length eL in addition to lens[seg],
     and take the part of each segment starting form offset, with nL
     subtracted.

     Returns the length of the shifted series.
  */
  int i=0, n=0, seg=0, cL=0;
  cL=0;
  for(seg=0; seg<n_segs; seg++) {
    n = lens[seg] + eL - nL;
    if(n > 0) {
      memcpy(outx+i, x+cL+offset, sizeof(ELEM)*n);
    } else {
      n = 0;
    }
    cL += lens[seg] + eL;
    i += n;
  }
  return i;
}

/* test whether x and y are significantly associated (e.g. by correlation),
   and return a score (larger means more significant) and the test_value
 */
typedef void (*test_func)(int L, ELEM x[], ELEM y[], double* score,double* test_value);

int test_link(int L, ELEM x[], ELEM y[], int min_delay, int max_delay, double st,
	      test_func f, int n_segs, int lens[],
	      links out_links[], ELEM buf[]) {
  /* x and y have length L, representing two time-series.
     st is the score threshold, below which the link is not significant.
     f is a function that tests whether two vectors have significant association.
     lens contains n_segs lengths of the segments in x and y.
     out_links should have room for up to max_delay links.
     buf is for temporary use, and should have length 2L.

     test whether x --> y.
     return the number of links added.
   */
  int i=0, sL=0, n_links=0;
  double score=0, test_value=0;

  /* debug */
  /*
  printf("==== test_link\n");
  printf("lens"); pr_ints(n_segs,lens);
  printf("x"); pr_doubles(L,x);
  printf("y"); pr_doubles(L,y);
  */

  for(i=min_delay; i<=max_delay; i++) {
    sL = shift_by(x,n_segs,lens, 0, i, buf);
    if(sL <= 5) continue; /* not long enough */
    shift_by(y,n_segs,lens, i, i, buf+L);

    /* debug */
    /*
    printf("** at delay %d\n", i);
    */
    /*
    printf("at delay %d, x", i); pr_doubles(sL,buf);
    printf("at delay %d, y", i); pr_doubles(sL,buf+L);
    */

    f(sL, buf,buf+L, &score,&test_value);

    if(score >= st) { /* add a link */
      out_links[n_links].delay = i;
      out_links[n_links].score = score;
      out_links[n_links].test_value = test_value;
      n_links++;
    }
  }
  return n_links;
}

int infer_grn1(array2d* d, double st, int min_delay, int max_delay, test_func f,
	       int n_segs, int lens[], char* forbidden,
	       int out_nL[], links* out_Ls[],
	       links buf_links[], ELEM buf[]) {
  /* Stage 1 of DCLINDE, infers the time lag and pairwise association.
     d is n by g array, the n_segs segments of expression data in column major format,
     where rows are time points, and columns are genes.
     The first lens[0] rows are segment 1,
     the next lens[1] rows are segment 2, and so on.
     st is the score threshold, below which the link is not significant.
     f is a function that tests whether two vectors have significant association.

     out_nL should have length g*g, which will be treated in row major
     format, so that out_nL[j+i*g] is the number of links of i --> j.
     out_Ls is similar to out_nL, but holds the pointer to the links,
     i.e. out_Ls[j+i*g] is the array of links of i --> j.

     if forbidden is non-NULL, then if forbidden[i*g + j] is 1, then the link i->j is forbidden.

     buf_links should have room for up to g*g*max_delay links,
     and the pointers in out_Ls will point to here.

     buf should have room for up to 2n.

     Returns the number of links added.
   */
  int n=0, g=0, i=0, j=0, k=0, nL=0, n_links=0;
  n = d->rows;
  g = d->cols;

  for(i=0; i<g; i++) {
    for(j=0; j<g; j++) {
      /* test i --> j */
      /* debug */
      /*
      printf("** test %d -> %d\n", i,j);
      */

      /* allow self loop, but only for positive delays */
      if((forbidden!=NULL && forbidden[i*g+j]==1)) {
	out_nL[j+i*g] = 0;
	continue;
      }
      nL = test_link(n, get_col(d,i),get_col(d,j),
		     i==j ? (min_delay>0 ? min_delay : 1) : min_delay,
		     max_delay, st,
		     f, n_segs,lens, buf_links+n_links,buf);
      out_nL[j+i*g] = nL;
      out_Ls[j+i*g] = buf_links+n_links;
      /* still need to fill in the from and to */
      for(k=0; k<nL; k++) {
	buf_links[n_links + k].from = i;
	buf_links[n_links + k].to = j;
      }
      n_links += nL;
    }
  }
  return n_links;
}

static int best_link(int n, links Ls[]) {
  double s = Ls[0].score;
  int i=1, idx=0;
  for(i=1; i<n; i++) {
    if(Ls[i].score > s) {
      s = Ls[i].score;
      idx = i;
    }
  }
  return idx;
}
int remove_dup_links(int g, int nL[], links* Ls[]) {
  /* both nL and Ls are g by g in row major format.
     nL[j+i*g] stores the number of links of i --> j.
     Ls[j+i*g] stores the array of links of i --> j, with different delays.
     modify Ls such that for each i --> j, there are at most one delay,
     by choosing the one with the highest score.

     Returns the number of links removed.
   */
  int i=0,j=0, k=0, n=0, n_removed=0;
  links* p=NULL;

  for(i=0; i<g; i++) {
    for(j=0; j<g; j++) {
      n = nL[j+i*g];
      if(n <= 1) continue;
      p = Ls[j+i*g];
      k = best_link(n, p);
      if(k > 0) p[0] = p[k];
      n_removed += n-1;
      nL[j+i*g] = 1;
    }
  }
  return n_removed;
}

int remove_zero_delay_links(int n_links, links* p_links[]) {
  /* p_links are the n_links pointers to links.
     modify p_links to remove those with zero delays.
     returns the number of links after removal.
   */
  int i=0, j=0;
  for(i=0; i<n_links; i++) {
    if(p_links[i]->delay == 0) continue;
    p_links[j] = p_links[i];
    j++;
  }
  return j;
}

int remove_dup_links2(int g, int nL[], links* Ls[], int n_links, links* p_links[]) {
  /* both nL and Ls are g by g in row major format.
     nL[j+i*g] stores the number of links of i --> j.
     Ls[j+i*g] stores the array of links of i --> j, with different delays.

     p_links are the n_links pointers to links.
     modify Ls such that for each i --> j, there are at most one delay,
     by choosing the one with the highest score.

     Returns the number of links after removal.
   */
  int i=0,j=0, k=0, n=0;
  int from=0, to=0;
  links* p=NULL;

  /* compact p_links, mark and remove dup links, adjust nL */
  for(i=0, j=0; i<n_links; i++) {

    /* debug */
    /*
    printf("compacting i: %d\tp_links[i]:%p\tfrom: %d\tto:%d\tdelay:%d\t%g\n", i,p_links[i],p_links[i]->from,p_links[i]->to,p_links[i]->delay, p_links[i]->score);
    */

    if(p_links[i]->delay < 0) { /* marked to be removed */
      p_links[i] = NULL;
      continue;
    }

    /* keep this link */
    from = p_links[i]->from;
    to = p_links[i]->to;
    n = nL[to + from*g];

    if(i!=j) p_links[j] = p_links[i];
    j++; /* to keep the link, increment no matter whether i==j */

    if(n > 1) { /* mark other delays to be removed */
      p = Ls[to + from*g];
      for(k=0; k<n; k++) {
	if(p+k != p_links[i])
	  p[k].delay = -1; /* mark it as not used */
      }
      nL[to + from*g] = 1;
    }
  }
  /* adjust Ls */

  /* debug */
  /*
  printf("number of links after no dup: %d\n", j);
  */

  for(i=0; i<j; i++) {
    from = p_links[i]->from;
    to = p_links[i]->to;

    /* debug */
    /*
    printf("i: %d\tp_links[i]:%p\tfrom: %d\tto:%d\tdelay:%d\n", i,p_links[i],from,to,p_links[i]->delay);
    printf("\tLs[to+from*g]: %p\tnL[to+from*g]:%d\n", Ls[to + from*g], nL[to + from*g]);
    fflush(stdout);
    */

    if(p_links[i] != Ls[to + from*g]) {
      Ls[to + from*g][0] = *(p_links[i]);
      p_links[i] = Ls[to + from*g];
    }
  }

  return j;
}

static void print_one_link(FILE* f, int from, int to, int delay) {
  fprintf(f, "To: %d From: %d Delay: %d\n", to, from, delay);
}

static void print_n_links(FILE* f, int n_links, links links[]) {
  int i=0;
  for(i=0; i<n_links; i++)
    print_one_link(f, links[i].from, links[i].to, links[i].delay);
}

void print_link(FILE* f, links* L) {
  print_one_link(f, L->from, L->to, L->delay);
}
void print_only_links(FILE* f, int n_links, links* p_links[]) {
  int i=0;
  for(i=0; i<n_links; i++) {
    print_link(f, p_links[i]);
  }
}
void print_links(FILE* f, int n_links, links* p_links[]) {
  fprintf(f, "number of links: %d\n", n_links);
  print_only_links(f, n_links, p_links);
}

/***************************************************************/
int shift_c_vector(array2d* d, int n_segs, int lens[],
		   int ix, int iy, int lag,
		   int nZ, int iZ[], int dZ[],
		   ELEM out_x[], ELEM out_y[], ELEM out_z[]) {
  /* To shift the x, y and Z by lags as appropriate.
     d is the data in column major format, n by g, with n_segs segments,
     each with lengths in lens.
     ix and iy are the 0-based column index into d.
     lag is the delay of x --> y.

     iZ contains nZ 0-based column indices into d,
     dZ contains the corresponding delays with respect to x.

     output the shifted parts into out_x and out_y for x and y.
     output the shifted parts into out_z for Z, column by column,
     so the total length of out_z will be length(out_x)*nZ.

     Assume the segments are long enough such that shifting does not
     make the segments crossing each other.

     Returns the length of the extracted out_x and out_y.
   */
  int md=0, k=0, L=0, sL=0;
  /* md = min(0, lag, dZ), each delay is to be subtracted from md */
  md = min_int_of(nZ, dZ);
  if(md > 0) md = 0;
  if(md > lag) md = lag;

  /* L = max(0 - md, lag - md, dZ - md).
     L is the length threshold, below which a segment will be shifted to empty
   */
  L = max_int_of(nZ,dZ) - md;
  if((0-md) > L) L = 0 - md;
  if((lag-md) > L) L = lag - md;

  /* take x, original delay is 0 */
  sL = shift_by(get_col(d, ix), n_segs,lens,  0 - md, L, out_x);
  if(sL <= 0) return 0;
  /* take y, original delay is lag */
  shift_by(get_col(d, iy), n_segs,lens,  lag - md, L, out_y);
  /* take Z, oringinal delays in dZ */
  for(k=0; k<nZ; k++) {
    shift_by(get_col(d, iZ[k]), n_segs,lens,  dZ[k] - md, L, out_z+k*sL);
  }
  return sL;
}

/* test whether x and y are significantly associated (e.g. by
   correlation), given Z, and return a score (larger means more
   significant) and the test_value.
   Z has nZ columns, each of length L.
 */
typedef void (*test_c_func)(int L, ELEM x[], ELEM y[], int nZ, ELEM Z[], double* score,double* test_value);

int test_c_link(array2d* d, int n_segs, int lens[],
		double st, test_c_func f,
		int ix, int iy, int lag,
		int n_neis, int neis[],
		int n_d_neis[], int* d_neis[],
		int n_z, int z_set[],
		int ibuf[], ELEM buf[]) {
  /* to do condition test of x --> y | Z.
     most parameters are same as in shift_c_vector().
     neis contains the n_neis indices of the neighbors of x --> y.
     n_d_neis contains the number of delays for each neighbor.
     d_neis contains the array of delays for each neighbor.
     z_set contains n_z indices into neis, z_set is the set of
     Z we are conditioning on.

     ibuf should have length 3g, where g=ncol(d),
     buf should have length (n_z+2)*nrow(d).

     st is score threshold, below which the test is insignificant,
     f is the testing function.

     Returns 0 if any of the test is insignificant.
     Returns 1 otherwise.
  */
  int i=0, n=0, g=0, L;
  int *iZ=NULL, *ic=NULL, *dZ=NULL;
  ELEM *x=NULL, *y=NULL, *Z=NULL;
  double score=0, test_value=0;

  n = d->rows;
  g = d->cols;

  iZ = ibuf;
  ic = iZ+g;
  dZ = ic+g;

  x = buf;
  y = x + n;
  Z = y + n;

  for(i=0; i<n_z; i++) /* become indices of columns into d */
    iZ[i] = neis[z_set[i]];

  /* debug */
  /*
  printf("test_c_link, ix:%d, iy:%d, iZ:", ix,iy); pr_ints(n_z,iZ); fflush(stdout);
  */


  /* try out the delays of Z */
  for(i=0; i<n_z; i++) ic[i] = 0; /* first set of delays, as counter */
  while(1) {
    /* get the delays */
    for(i=0; i<n_z; i++) dZ[i] = d_neis[z_set[i]][ic[i]];
    /* debug */
    /*
    printf("test_c_link, dZ:"); pr_ints(n_z,dZ); fflush(stdout);
    */

    L = shift_c_vector(d,n_segs,lens, ix,iy,lag, n_z,iZ,dZ, x,y,Z);
    if(L > 5) { /* ignore if too short */
      score = 0;
      f(L, x,y, n_z,Z, &score,&test_value);

      /* debug 
      printf("test_c_link, score: %f\n", score); fflush(stdout);
      */

      if(score < st) {
      /* test fail for whatever reason, or score not significant */
	return 0;
      }
    }
    /* increment the counter of the delays */
    for(i=n_z-1; i>=0; i--) {
      if(++ic[i] < n_d_neis[z_set[i]]) break; /* incremented, no overflow */
      /* one overflow, reset it, increment the next digit */
      ic[i] = 0;
    }
    if(i < 0) break; /* tried through the delays */
  }

  /* seems good */
  return 1;
}

/***************************************************************/
static int get_neighbors(int g, int ix, int iy, int nL[],
			 int* out_n_x, int* out_n_y,
			 int out_x[], int out_y[]) {
  /* neighbors of x and y except x and y themselves,
     nL[j+i*g] is the number of i --> j links.
     output the indices of the neighbors in out_x and out_y respectively,
     first the parents, then children, also ouput the counts.

     return the total number of neighbors.
   */
  int i=0, n_x=0,n_y=0;

  for(i=0; i<g; i++) {
    if(i==ix || i==iy) continue;
    if(nL[ix + i*g] + nL[i + ix*g] > 0) out_x[n_x++] = i;
    if(nL[iy + i*g] + nL[i + iy*g] > 0) out_y[n_y++] = i;
  }
  /* done */
  if(out_n_x) *out_n_x = n_x;
  if(out_n_y) *out_n_y = n_y;

  return n_x + n_y;
}

static int get_common_neighbors(int g, int ix, int iy, int nL[],
				int out[]) {
  /* common neighbors of x and y except x and y themselves,
     nL[j+i*g] is the number of i --> j links.
     output the indices of the common neighbors in out

     return the total number of common neighbors.
   */
  int i=0, n=0;

  /* into */
  for(i=0; i<g; i++) {
    if(i==ix || i==iy) continue;
    /* (i --> ix or ix --> i) && (i --> iy or iy --> i) */
    if((nL[ix + i*g] + nL[i + ix*g] > 0)
       && (nL[iy + i*g] + nL[i + iy*g] > 0))
      out[n++] = i;
  }

  return n;
}

static void get_the_delays_to_x(int g, int ix, int iy, int lag,
				int nL[], links* Ls[],
				int include_x,int include_y,
				int nZ, int iZ[],
				int out_ndZ[], int* out_dZ[], int ibuf[]) {
  /* get the delays of Z relative to x.
     g, ix, iy and nL have the same meaning as in get_common_neighbors().
     lag is the delay of x --> y.
     iZ contains nZ indices of the Z's.

     include_x indicates whether to include the delays of neighbors of x.
     include_y is similar, but for y.

     output the number of delays for each z in out_ndZ,
     output the delays for each z in out_dZ, and setup the pointer that
     point into ibuf.
     ibuf should be large enough to hold the needed delays.
   */
  int i=0, j=0, k=0, z=0, n=0;
  links* p=NULL;

  for(i=0; i<nZ; i++) {
    z = iZ[i];

    k = 0;
    if(include_x) {
      /* x --> z */
      n = nL[z + ix*g];
      p = Ls[z + ix*g];
      for(j=0; j<n; j++) ibuf[k++] = p[j].delay;
      /* z --> x */
      n = nL[ix + z*g];
      p = Ls[ix + z*g];
      for(j=0; j<n; j++) ibuf[k++] = -(p[j].delay);
    }

    if(include_y) {
      /* y --> z */
      n = nL[z + iy*g];
      p = Ls[z + iy*g];
      for(j=0; j<n; j++) ibuf[k++] = lag + (p[j].delay);
      /* z --> y */
      n = nL[iy + z*g];
      p = Ls[iy + z*g];
      for(j=0; j<n; j++) ibuf[k++] = lag - (p[j].delay);
    }
    /***/
    out_ndZ[i] = k;
    out_dZ[i] = ibuf;
    ibuf += k;
  }
}

static int enumerate_next_subset(int n, int size, int set[]) {
  /* set contains size numbers in 0 .. n-1, in strictly increasing order,
     representing a subset of size numbers in 0 .. n-1.
     Try to get to the next subset by increasing the numbers.
     Return 1 if OK. Return 0 if cannot enumerate the next subset anymore.
     set is modified in place
   */
  int i=0, j=0, k=0, idx;
  for(i=0; i<size; i++) { /* to start from the end */
    idx = size-1-i;
    if(++set[idx] < n-i) { /* the last one can hold a max of n-1, second last can hold a max of n-2 */
      /* OK, reinitialize the rest */
      k = set[idx];
      for(j=idx+1; j<size; j++) set[j] = ++k;
      return 1;
    }
  }
  /* tried through */
  return 0;
}

int filter_links(array2d* d, int n_segs, int lens[],
		 int nL[], links* Ls[],
		 int n_links, links* p_links[],
		 double st, int max_n, test_c_func f, int pruning) {
  /* Stage 2 of DCLINDE.
     d is n by g array, the n_segs segments of expression data in column major format,
     where rows are time points, and columns are genes.
     The first lens[0] rows are segment 1,
     the next lens[1] rows are segment 2, and so on.

     st is the score threshold, below which a link is removed on conditional test.
     max_n is the maximum number of neighbors to condition on.
     f is the function for doing conditional test.

     nL[j+i*g] is the number of i --> j links,
     Ls[j+i*g] is the pointer to array of i --> j links.
     p_links contains the pointers to the n_links.

     pruning is either PRUNING_ALL (choose from all neighbors) or
     PRUNING_COMMON (choose only from common neighbors) for condition
     on neighbors.

     nL, Ls and p_links will be modified by removing links.
     Returns the number of links after removal.
   */
  int cn=0, i=0, j=0, is_end=0, n=0,g=0;
  int from=0, to=0, lag=0;
  int n_neis=0, n_x=0,n_y=0, include_x=0,include_y=0;
  int *neis=NULL, *nei_x=NULL,*nei_y=NULL,
    *idx=NULL, *iZ=NULL, *ibuf=NULL,
    *n_d_neis=NULL, *d_neis_buf=NULL;
  int** d_neis = NULL;
  ELEM *buf=NULL;

  if(n_links <= 0) return 0; /* no links, no need to do anything */

  n = d->rows;
  g = d->cols;

  i = 8*g + n_links;

  neis = Malloc(i, int, "filter_links()");
  memset(neis, 0, sizeof(int)*i);
  nei_x = neis;            /* g */
  nei_y = nei_x+g;         /* g */
  idx = nei_y+g;           /* g */
  iZ = idx+g;              /* g */
  ibuf = iZ+g;             /* 3g */
  n_d_neis = ibuf+3*g;       /* g */
  d_neis_buf = n_d_neis+g; /* len n_links */

  d_neis = Malloc(g, int*, "filter_links()");
  memset(d_neis, 0, sizeof(int*)*g);

  buf = Malloc((max_n+2)*n, ELEM, "filter_links()");
  memset(buf, 0, sizeof(ELEM)*(max_n+2)*n);

  /* condition on more and more neighbors */
  for(cn=1; !is_end && cn<=max_n; cn++) {
    is_end = 1;

    /* debug */
    /*
    printf("====== iteration %d\n", cn);
    */

    /* n_links might be updated after the following loop */
    for(i=n_links-1; i>=0; i--) {
      if(p_links[i] == NULL) continue;

      /* debug */
      /*
      printf("p_link[%d] (%p) before adjust: to: %d from: %d lag: %d\n", i, p_links[i], p_links[i]->to, p_links[i]->from, p_links[i]->delay);
      */

      /* test links with lower scores first */
      from = p_links[i]->from;
      while(from < 0) { /* marked to update p_links */
	/* may need to update a few times */
	p_links[i] += from;
	from = p_links[i]->from;
      }
      to = p_links[i]->to;
      lag = p_links[i]->delay;

      /* debug */
      /*
      printf("p_link[%d] (%p) after adjust: to: %d from: %d lag: %d\n", i, p_links[i], p_links[i]->to, p_links[i]->from, p_links[i]->delay);
      printf("== p_links[%d]: to: %d from: %d lag: %d\n", i, to,from,lag);
      */

      /* determine neighbors to condition on */
      include_x=1;
      include_y=1;

      if(pruning == PRUNING_COMMON) {
	n_neis = get_common_neighbors(g, from,to, nL, neis);
      } else {
	get_neighbors(g, from,to, nL, &n_x,&n_y, nei_x,nei_y);
	/* debug */
	/*
	printf("= prune all\n");
	printf("nei_x"); pr_ints(n_x,nei_x);
	printf("nei_y"); pr_ints(n_y,nei_y);
	*/

	if(n_x < n_y) {
	  n_neis = n_x;
	  include_y = 0;
	} else {
	  n_neis = n_y;
	  memcpy(neis, nei_y, n_neis*sizeof(int));
	  include_x = 0;
	}
      }

      /* debug */
      /*
      printf("== neis");
      pr_ints(n_neis, neis);
      */

      if(n_neis < cn) continue;
      /* need to check */
      is_end = 0;

      get_the_delays_to_x(g, from,to,lag, nL,Ls, include_x,include_y,
			  n_neis,neis,
			  n_d_neis,d_neis, d_neis_buf);

      /* debug */
      /*
      printf("n_d_neis"); pr_ints(n_neis, n_d_neis);
      printf("d_neis"); pr_int_arr(n_neis, n_d_neis,d_neis);
      */

      /* generate combinations of neighbors */
      for(j=0; j<cn; j++) idx[j] = j; /* first subset, 0-based index */
      do {
	/* debug */
	/*
	printf("subset"); pr_ints(cn,idx);
	*/

	/* do conditional test */
	if(!test_c_link(d,n_segs,lens, st,f, from,to,lag,
			n_neis,neis, n_d_neis,d_neis,
			cn,idx, ibuf,buf)) {
	  /* remove it */
	  {
	    /* debug */
	    /*
	    printf("=== remove p_links[%d]\n", i);
	    */

	    int m = nL[to + from*g];
	    links* pL = Ls[to + from*g] + m-1;

	    /* debug */
	    /*
	    printf("pL (%p): to: %d from: %d lag: %d\n", pL, pL->to, pL->from, pL->delay);
	    printf("pL - p_links[i]: %d\n", (pL - p_links[i]));
	    */

	    if(pL != p_links[i]) {
	      *(p_links[i]) = *pL;
	      /* to mark the last link to update p_links when encountered */
	      pL->from = -(pL - p_links[i]);
	    }

	    /* debug */
	    /*
	    printf("p_link[%d] (%p) after: to: %d from: %d lag: %d\n", i, p_links[i], p_links[i]->to, p_links[i]->from, p_links[i]->delay);
	    printf("pL (%p) after: to: %d from: %d lag: %d\n", pL, pL->to, pL->from, pL->delay);
	    */

	    nL[to + from*g]--;
	    p_links[i] = NULL;
	  }
	  break;
	}
      } while(enumerate_next_subset(n_neis, cn, idx));
    }
    /* compact p_links if needed */
    for(i=0, j=0; i<n_links; i++) {
      if(p_links[i] != NULL) {
	/* correct the possibly -1 from */
	while(p_links[i]->from < 0) p_links[i] += p_links[i]->from;

	if(i != j) p_links[j] = p_links[i];
	j++; /* to keep the link, increment no matter whether i==j */
      }
    }
    n_links = j;

    /* debug */
    /*
    printf("=== links after iteration %d\n", cn);
    print_links(stdout, n_links, p_links);
    */
  }
  /* clean up */
  free(neis);
  free(d_neis);
  free(buf);
  return n_links;
}

/***************************************************************/
int decreasing_link_score(const void* a, const void* b) {
  links *A,*B;
  A = *(links**)a;
  B = *(links**)b;

  if(A->score > B->score) return -1;
  if(A->score < B->score) return 1;

  return A->delay - B->delay;
}
static int sort_links(int g, int nL[], links* Ls[], links* out_links[]) {
  /* nL[j+i*g] is the number of i --> j links,
     Ls[j+i*g] is the pointer to array of i --> j links.
     sort the links by decreasing score and put the pointers in out_links.

     Returns the number of links.
   */
  int i=0, j=0, k=0, n=0, n_links=0;
  links* p=NULL;

  for(i=0; i<g; i++) {
    for(j=0; j<g; j++) {
      n = nL[j+i*g];
      if(n <= 0) continue;
      p = Ls[j+i*g];
      for(k=0; k<n; k++) {
	out_links[n_links++] = p+k;
      }
    }
  }
  /* sort it */
  qsort(out_links, n_links, sizeof(links*), decreasing_link_score);
  /**/
  return n_links;
}

int dclinde(array2d* d, int n_segs, int lens[],
	    double st1, double st2, int min_delay, int max_delay, int max_n,
	    int pruning, int is_one_delay, int is_no_dup, int is_keep_zero_delay,
	    char* forbidden, links** out_links,
	    char* out_grn, char* out_s1_grn) {
  /* main function of DCLINDE.
     d is n by g array, the n_segs segments of expression data in column major format,
     where rows are time points, and columns are genes.
     The first lens[0] rows are segment 1,
     the next lens[1] rows are segment 2, and so on.

     st1 and st2 are the score thresholds for stage 1 and 2 respectively.
     max_n is the maximum number of neighbors to condition on.

     method is either METHOD_PCOR (partial correlation) or METHOD_MI (mutual information).
     pruning is either PRUNING_ALL (choose from all neighbors) or PRUNING_COMMON (choose only from common neighbors) for condition on neighbors in stage 2.

     if is_one_delay is true, remove duplicate links after stage 1.
     if is_keep_zero_delay is true, keep links with zero delay after stage 2 and before removing dups.
     if is_no_dup is true, remove duplicate links after stage 2.

     if forbidden is non-NULL, then if forbidden[i*g + j] is 1, then the link i->j is forbidden for 0 <= i,j < ng.
     if out_links is non-NULL and there are links, *out_links will be written with a pointer to the final links, and the buffer should be freed after use.

     if out_grn is non-NULL, output the final GRN to file out_grn.
     if out_s1_grn is non-NULL, output the stage 1 GRN to file out_s1_grn.

     returns the number of final links.
   */
  int n_links=0, g=0,n=0;
  int* nL=NULL;
  links** Ls=NULL;
  links** p_links=NULL;
  links* buf_links=NULL;
  ELEM* buf=NULL;
  FILE* tmpf=NULL;

  n = d->rows;
  g = d->cols;

  /* prepare buffers */
  nL = Malloc(g*g, int, "infer_grn()");
  memset(nL, 0, sizeof(int)*g*g);

  Ls = Malloc(g*g, links*, "infer_grn()");
  memset(Ls, 0, sizeof(links*)*g*g);

  p_links = Malloc(g*g*(1+max_delay), links*, "infer_grn()");
  memset(p_links, 0, sizeof(links*)*g*g*(1+max_delay));

  buf_links = Malloc(g*g*(1+max_delay), links, "infer_grn()");
  memset(buf_links, 0, sizeof(links)*g*g*(1+max_delay));

  buf = Malloc(2*n, ELEM, "infer_grn()");
  memset(buf, 0, sizeof(ELEM)*2*n);

  /* stage 1 */
  report_stop_watch();
  printf("==== Stage 1\n");
  infer_grn1(d, st1, min_delay, max_delay,
	     test_func_gtest,
	     n_segs, lens, forbidden, nL,Ls, buf_links, buf);
  if(is_one_delay) remove_dup_links(g, nL, Ls);
  n_links = sort_links(g, nL, Ls, p_links);

  printf("==== Initial GRN after stage 1:\n");
  print_links(stdout, n_links, p_links);
  printf("==== End stage 1 ==============\n");
  /* to out_s1_grn if applicable */
  if(out_s1_grn != NULL) {
    printf("=== output stage 1 GRN to file %s.\n", out_s1_grn);
    tmpf = fopen(out_s1_grn, "wt");
    if(tmpf == NULL) {
      fprintf(stderr, "Error: Cannot open file %s to output stage 1 GRN.\n",
	      out_s1_grn);
    } else {
      print_only_links(tmpf, n_links, p_links);
      fclose(tmpf);
      tmpf = NULL;
    }
  }

  report_stop_watch();
  fflush(stdout);
  /* stage 2 */
  printf("==== Stage 2\n");
  n_links = filter_links(d, n_segs, lens, nL, Ls, n_links, p_links,
			 st2, max_n,
			 test_c_func_gtest,
			 pruning);
  /***/
  if(!is_keep_zero_delay) n_links = remove_zero_delay_links(n_links, p_links);
  /* note that nL and Ls may have more links than in p_links,
     so the main reference are n_links and p_links.
   */
  if(is_no_dup) n_links = remove_dup_links2(g, nL, Ls, n_links, p_links);

  printf("==== \n");
  print_links(stdout, n_links, p_links);
  printf("==== End stage 2 ============\n");
  /* to out_grn if applicable */
  if(out_grn != NULL) {
    printf("=== output GRN to file %s.\n", out_grn);
    tmpf = fopen(out_grn, "wt");
    if(tmpf == NULL) {
      fprintf(stderr, "Error: Cannot open file %s to output GRN.\n",
	      out_grn);
    } else {
      print_only_links(tmpf, n_links, p_links);
      fclose(tmpf);
      tmpf = NULL;
    }
  }
  fflush(stdout);

  /* copy the links to out_links if applicable */
  if(out_links != NULL && n_links>0) {
    {
      links* outL = NULL;
      int i=0;
      outL = Malloc(n_links, links, "dclinde()");
      for(i=0; i<n_links; i++) {
	outL[i] = *(p_links[i]);
      }
      *out_links = outL;
    }
  }

  /* clean up */
  free(nL);
  free(Ls);
  free(p_links);
  free(buf_links);
  free(buf);

  return n_links;
}

typedef int (*pf_globalMIT_wrapper)(array2d* input_data, int n, int Dim, int n_segs, int lens[],
				    int order, int i_start, int i_end, double Alpha, int allow_self_link,
				    char* forbidden, links** out_links);
int globalMIT_exe_wrapper(array2d* input_data, int n, int Dim, int n_segs, int lens[],
			  int order, int i_start, int i_end, double Alpha, int allow_self_link,
			  char* forbidden, links** out_links); /* defined elsewhere */
int globalMIT_fixedOrder_exe_wrapper(array2d* input_data, int n, int Dim, int n_segs, int lens[],
				     int order, int i_start, int i_end, double Alpha, int allow_self_link,
				     char* forbidden, links** out_links); /* defined elsewhere */

int globalMIT_exe(array2d* d, int n_segs, int lens[],
		  int max_delay, int i_start, int i_end, double alpha, int allow_self_link,
		  pf_globalMIT_wrapper wrapper_func,
		  char* forbidden, links** out_links, char* out_grn) {
  /*
      if forbidden is non-NULL, then if forbidden[i*dim + j] is 1, then the link i->j is forbidden for 0 <= i,j < dim.
      where dim is d->cols.
  */
  int n_links=0, i=0;
  links* Ls=NULL;
  FILE* tmpf=NULL;

  n_links = wrapper_func(d, d->rows, d->cols, n_segs, lens,
			 max_delay, i_start, i_end < 0 ? (d->cols)-1 : i_end,
			 alpha, allow_self_link, forbidden, &Ls);

  /* to out_grn if applicable */
  if(out_grn != NULL) {
    printf("=== output GRN to file %s.\n", out_grn);
    tmpf = fopen(out_grn, "wt");
    if(tmpf == NULL) {
      fprintf(stderr, "Error: Cannot open file %s to output GRN.\n",
	      out_grn);
    } else {
      for(i=0; i<n_links; i++)
	print_link(tmpf, Ls+i);
      fclose(tmpf);
      tmpf = NULL;
    }
  }
  fflush(stdout);

  if(out_links) {
    *out_links = Ls;
  } else {
    free(Ls);
  }

  return n_links;
}

/***************************************************************/
static int int_prod_of(int n, int x[]) {
  int r=1, i=0;
  for(i=0; i<n; i++)
    r *= x[i];
  return r;
}

static int states_to_index(int L, int states[], int n_states[]) {
  /* index of the states vector, with respect to the number of states */
  int i=0, r=1, s=0;
  for(i=0; i<L; i++) {
    s += states[i]*r;
    r *= n_states[i];
  }
  return s;
}

static void index_to_states(int idx, int L, int n_states[], int out_states[]) {
  /* convert the index idx back to states in out_states */
  int i=0;
  for(i=0; i<L; i++) {
    out_states[i] = idx % n_states[i];
    idx /= n_states[i];
  }
}

static void print_index_as_states(FILE* f, int idx, int L, int n_states[]) {
  /* same as index_to_states(), but print it directly */
  int i=0;
  for(i=0; i<L; i++) {
    fprintf(f," %d", idx % n_states[i]);
    idx /= n_states[i];
  }
}

static void normalize_dist(int L, double x[]) {
  int i=0;
  double s=0;

  for(i=0; i<L; i++) s += x[i];
  if(s < 1e-12) s=1e-12;

  for(i=0; i<L; i++) x[i] /= s;
}

typedef struct s_node {
  int n_states;
  int n_parents;
  /* one single buffer parents, delays and parent_states, so free parents only */
  int* parents; /* own this, length n_parents, 0-based indices of parents */
  int* delays;  /* own this, length n_parents, delays of each parent link */
  int* parent_n_states; /* own this, length n_parents, number of states for each parent */
  int n_cpt;    /* number of conditional distributions, should be same as int_prod_of(n_parents, parent_n_states) */
  /* need to free cpt */
  double* cpt; /* own this, every n_states is a conditional
		  distribution, total int_prod_of(n_parents,
		  parent_states) conditional distributions. */
} node;

void destroy_node(node* n) {
  if(n==NULL) return;
  if(n->parents) free(n->parents); /* this also frees delays and parent_states */
  if(n->cpt) free(n->cpt);
  memset(n, 0, sizeof(node));
}

void free_nodes(int n_nodes, node* ns) {
  /* also frees ns itself */
  int i=0;
  for(i=0; i<n_nodes; i++) destroy_node(ns+i);
  if(ns) free(ns);
}

static void pr_ints_to_f(FILE* f, int n, int x[]) {
  int i=0;
  for(i=0; i<n; i++)
    printf(" %d", x[i]);
}

static void pr_doubles_to_f(FILE* f, int n, double x[]) {
  int i=0;
  for(i=0; i<n; i++)
    printf(" %g", x[i]);
}

static void print_nodes(FILE* f, int n_nodes, node* ns) {
  /* all node indices are 0-based */
  int i=0, j=0, n_ps=0;
  fprintf(f, "Number of nodes: %d\n", n_nodes);
  for(i=0; i<n_nodes; i++) {
    fprintf(f, "\nNode %d:\n", i);
    fprintf(f, "n.states: %d\n", ns[i].n_states);
    fprintf(f, "n.parents: %d\n", ns[i].n_parents);
    fprintf(f, "parents:");
    pr_ints_to_f(f, ns[i].n_parents, ns[i].parents);
    fprintf(f, "\n");
    fprintf(f, "delays:");
    pr_ints_to_f(f, ns[i].n_parents, ns[i].delays);
    fprintf(f, "\n");
    fprintf(f, "parents n.states:");
    pr_ints_to_f(f, ns[i].n_parents, ns[i].parent_n_states);
    fprintf(f, "\n");

    fprintf(f, "CPT:\n");
    n_ps = int_prod_of(ns[i].n_parents, ns[i].parent_n_states);
    for(j=0; j<n_ps; j++) {
      print_index_as_states(f, j, ns[i].n_parents, ns[i].parent_n_states);
      fprintf(f, " : ");
      pr_doubles_to_f(f, ns[i].n_states, ns[i].cpt + j*(ns[i].n_states));
      fprintf(f, "\n");
    }

    fprintf(f, "End Node %d\n", i);
  }
}

static node* links_to_nodes(int ng, int n_links, links* links,
			    array2d* data, int n_data, int lens[]) {
  /* allocates and fill in an array of ng nodes, with info from the links and data.

     data is used for estimating the CPT. data is assumed in column
     major format, with ng columns.
   */
  int i=0, j=0, k=0, m=0, from=0, to=0;
  int seg=0, sL=0, md=0, n_ps=0, *ps=NULL;
  node* ns=NULL;
  ns = Malloc(ng, node, "links_to_nodes()");
  memset(ns, 0, sizeof(node)*ng);

  m = data->rows;

  /* the number of states for each node */
  for(i=0; i<ng; i++)
    ns[i].n_states = 1+max_int_of(m, get_col(data,i));

  /* number of parents */
  for(i=0; i<n_links; i++) {
    to = links[i].to;
    ns[to].n_parents++;
  }

  /* allocate some buffers */
  for(i=0; i<ng; i++) {
    /* parents holds the buffer to parents, delays and parent_n_states,
       each of length n_parents. */
    ns[i].parents = Malloc(3*ns[i].n_parents, int, "links_to_nodes()");
    ns[i].delays = ns[i].parents + ns[i].n_parents;
    ns[i].parent_n_states = ns[i].delays + ns[i].n_parents;

    /* re-use n_parents as counter below */
    ns[i].n_parents = 0;
  }

  /* parents, delays, and parent_n_states */
  for(i=0; i<n_links; i++) {
    from = links[i].from;
    to = links[i].to;
    j = ns[to].n_parents;
    ns[to].parents[j] = from;
    ns[to].delays[j] = links[i].delay;
    ns[to].parent_n_states[j] = ns[from].n_states;

    ns[to].n_parents++;
  }

  /* CPT */
  for(i=0; i<ng; i++) {
    n_ps = int_prod_of(ns[i].n_parents, ns[i].parent_n_states);
    ns[i].n_cpt = n_ps;
    ns[i].cpt = Malloc(n_ps*(ns[i].n_states), double, "links_to_nodes()");
    for(j=n_ps*(ns[i].n_states)-1; j>=0; j--) {
      /* pseudo-count, so that zero counts become uniform distribution */
      ns[i].cpt[j] = 1e-6;
    }
    /* accumulate counts */
    ps = reserve_temp_buffer(ns[i].n_parents); /* no need to release here */
    md = ns[i].n_parents > 0 ? max_int_of(ns[i].n_parents, ns[i].delays) : 0;
    /* segment by segment */
    for(seg=0, sL=0; seg<n_data; seg++) {
      for(j=md; j<lens[seg]; j++) {
	/* parent states */
	for(k=0; k<ns[i].n_parents; k++) {
	  ps[k] = caref(data, sL+j-ns[i].delays[k], ns[i].parents[k]);
	}
	k = states_to_index(ns[i].n_parents, ps, ns[i].parent_n_states);
	/* count */
	ns[i].cpt[k*ns[i].n_states + caref(data,sL+j,i)] += 1;
      }
      sL += lens[seg];
    }
    /* normalize */
    for(j=0, k=0; j<n_ps; j++) {
      normalize_dist(ns[i].n_states, ns[i].cpt+k);
      k += ns[i].n_states;
    }
  }

  return ns;
}

static void calc_and_print_CPT(int ng, int n_links, links* links,
			       array2d* data, int n_data, int lens[]) {
  node* ns=NULL;
  ns = links_to_nodes(ng, n_links, links, data, n_data,lens);
  print_nodes(stdout, ng, ns);

  free_nodes(ng, ns);
}

int decreasing_double(const void* a, const void* b) {
  double ax=0, ay=0;
  ax = *(double*)a;
  ay = *(double*)b;
  if(ax > ay) return -1;
  if(ax < ay) return 1;
  return 0;
}

static double estimate_bias(int n_states, int n_distributions, double* cpt, double* buf) {
  /* cpt has n_distributions, each of length n_states.
     estimate the bias in an ad-hoc way.
     buf is for temporary use, should have length n_distributions.
   */
  int i=0;
  double* p=NULL;

  p = cpt;
  for(i=0; i<n_distributions; i++, p+=n_states) {
    buf[i] = max_of(n_states, p);
  }

  /* use the median of max of the conditional distributions */
  qsort(buf, n_distributions, sizeof(double), decreasing_double);
  return buf[n_distributions/2];
}

/***************************************************************/
/* A simple heuristic sequential clustering */
/* for recording the series included in each cluster */
typedef struct s_series {
  int id;
  int idx; /* column index */
  int offset;
  double variance; /* variance of the column */
  struct s_series* next;
} series;

typedef struct s_cluster {
  int id;
  int n_series; /* number of series merged to here */
  int capacity; /* length of buffer of mean and n, divided into segments */
  int s,e; /* 0-based start (inclusive) and end index (exclusive) of the section used */
  int* val; /* the value of each time point, own this */
  series* merged; /* not own this */
  struct s_cluster* next;
} cluster;

cluster* new_cluster(int id, int n_series, int capacity, int L,
		     int* init_val, int init_offset,
		     series* se,
		     cluster* next) {
  /* allocate a new cluster with the provided initial information (deep copy).
     If init_val is NULL, initialize to 0, and init_n is ignored and n set to 0.
     If init_val is non-NULL, copy its content as the mean; and if init_n is
     provided, it is copied as n, otherwise 1 is used.
  */
  cluster* c = NULL;

  c = Malloc(1,cluster, "new_cluster()");
  c->id = id;
  c->n_series = n_series;
  c->capacity = capacity;
  c->val = Malloc(capacity, int, "new_cluster()");
  memset(c->val, 0, sizeof(int)*capacity);

  assert(L <= capacity);
  c->s = init_offset;
  c->e = c->s + L;

  if(init_val != NULL) {
    memcpy(c->val + c->s, init_val, sizeof(int)*L);
  }

  c->merged = se; /* assume se is properly null-terminated */
  c->next = next;
  return c;
}

cluster* single_cluster(int id, int capacity, int L, int* init_val, int init_offset, series* se, cluster* next) {
  cluster* r=NULL;
  r = new_cluster(id, 1, capacity, L, init_val, init_offset, se, next);
  /* se->id is assumed to have been set properly */
  se->offset = r->s;
  se->next = NULL;
  return r;
}

void free_clusters(cluster* cs) {
  cluster* p=NULL;
  while(cs != NULL) {
    p = cs;
    cs = cs->next;
    if(p->val) free(p->val);
    free(p);
  }
}

/****************************************************************************/
double association_to_cluster(cluster* c, int n_segs, int lens[],
			      int L, int es[], int* offset, int buf[]) {
  /* calculate the max. association of es (of length m) to cluster c,
     return the association, and output the corresponding delay.
     es is of length L, which should be same as sum of lens.
     Assume es is wholly witin the capacity of c.
     The offset is relative to c.

     buf is for temporary use, should have length 2*L.
   */
  int s=0, best_offset=0, dx=0,dy=0, nL=0;
  int md=0, sL=0;
  double best_assoc=0, tmp_c=0;

  assert(c != NULL);

  /* assume that the initial space c->s is the allowed maximum delay. */
  md = c->s;

  for(s=0; s<=2*md; s++) {
    /* consider the segments, which are packed together starting at c->s */
    if(s < (c->s)) {
      dx = 0;
      dy = c->s - s;
      nL = dy;
    } else {
      dx = s - c->s;
      dy = 0;
      nL = dx;
    }

    sL = shift_by(c->val + c->s, n_segs,lens, dx, nL, buf);
    shift_by(es, n_segs,lens, dy, nL, buf+L);
    /* use G2 test, and -log_10(p-value) as association score */
    test_func_gtest(sL, buf, buf+L, &tmp_c,NULL);

    if(tmp_c > best_assoc) {
      best_assoc = tmp_c;
      best_offset = s;
    }
  }

  if(offset) *offset = best_offset;
  return best_assoc;
}

cluster* closest_cluster(cluster* cs, int n_segs, int lens[],
			 int L, int es[], int buf[],
			 int* out_offset, double* out_assoc) {
  /*
    cs is a list of clusters, go through each one, get the one with
    highest association.  return the cluster if found, return NULL if
    none satisfy the criteria.

    es is of length L, which should be same as sum of lens.
    buf is for temporary use, should have length 2*L.

    Also output the corresponding offset to *out_offset and
    association to *out_assoc.
   */
  cluster* r=NULL;
  double best_assoc=0, tmp_c=0;
  int best_offset=0, tmp_d=0;

  assert(cs != NULL);

  r = NULL;
  for(; cs!=NULL; cs=cs->next) {
    tmp_c = association_to_cluster(cs, n_segs,lens, L, es, &tmp_d, buf);
    printf("Association to cluster %d is %f, offset is %d\n",
	   cs->id, tmp_c, tmp_d);
    if(tmp_c > best_assoc) {
      r = cs;
      best_assoc = tmp_c;
      best_offset = tmp_d;
    }
  }

  if(out_offset) *out_offset = best_offset;
  if(out_assoc) *out_assoc = best_assoc;

  return r;
}

cluster* add_one_series(cluster* cs, series* se, int n_segs, int lens[],
			int L, int es[], int max_delay, double threshold,
			int buf[]) {
  /* add es (of n_segs segments, with lengths in lens) to the list of
     clusters in cs, either to an existing and closest cluster, or
     being a new cluster itself.

     es is of length L, which should be same as sum of lens.
     buf is for temporary use, should have length 2*L.

     threshold is for testing if the association score is large enough to include to a cluster.
     If the series is not close enough, es will become a new cluster.
     Return the modified linked list.
   */
  cluster* p=NULL;
  double assoc=0;
  int offset=0;

  /* at this stage, the segments are packed together with an offset.
     But the buffer of each cluster has room for each segment for shifting.
   */
  if(cs == NULL) {
    printf("New cluster 0\n");
    return single_cluster(0, L+n_segs*2*max_delay, L, es, max_delay, se, NULL);
  }

  p = closest_cluster(cs, n_segs,lens, L,es, buf, &offset,&assoc);

  if(p != NULL) {
    printf("Closest to cluster %d, association score: %f, offset: %d\n", p->id, assoc, offset);
  }

  if(p!=NULL && assoc >= threshold) { /* existing one */
    printf("Add to cluster %d\n", p->id);
    /* add only the series, but the val remains the first added series */
    p->n_series++;
    se->offset = offset;
    se->next = p->merged;
    p->merged = se;
    return cs;
  } else { /* new one */
    printf("New cluster %d\n", cs->id+1);
    return single_cluster(cs->id+1, L+n_segs*2*max_delay,L, es, max_delay, se, cs);
  }
}

cluster* find_clusters(array2d* d, int n_segs, int lens[],
		       int candidates[], int n_candidates,
		       int max_delay, double threshold, series s_buf[]) {
  /* assume d is in column major form, with n_segs segments, and the
     lengths of the segments are in lens.

     Consider the first n_candidates column indices in candidates for clustering.
     threshold is for testing if the association score is large enough to include to a cluster.

     s_buf is a buffer of series for use here, and will be referred to
     in the clusters.

     Return a list of clusters that have at least two series.
   */
  cluster *cs=NULL, *p=NULL, *r_cs=NULL;
  int m=0, i=0;
  int* buf=NULL;
  m = d->rows;

  for(i=0; i<n_candidates; i++) {
    s_buf[i].id = i;
    s_buf[i].idx = candidates[i];
    s_buf[i].offset = 0;
    s_buf[i].variance = 0; /* may need to change */
    s_buf[i].next = NULL;
  }

  buf = Malloc(2*m, int, "find_clusters()");
  /* sequentially go through the columns to find clusters */
  for(i=0; i<n_candidates; i++) {
    printf("\n** Consider series %d, column %d\n", i, candidates[i]);
    cs = add_one_series(cs, s_buf+i, n_segs,lens,
			m, get_col(d,candidates[i]),
			max_delay, threshold, buf);
  }
  /* print the clusters */
  printf("===== Resulting Clusters =====\n");
  for(p=cs; p!=NULL; p=p->next) {
    printf("Cluster %d: length %d, %d series\n", p->id, p->e - p->s, p->n_series);
    /* debug */
    /*
    for(i=0; i<(p->e); i++) {
      printf("\t%d: %d\n", i, p->val[i]);
    }
    printf("\n");
    */
  }

  /* filter out those with only one series */
  r_cs=NULL;
  while(cs != NULL) {
    /* debug */
    /*
    printf("cs: %x\tn_series: %d\tcapacity: %d\ts: %d, e: %d\tnext: %x\n",
	   cs, cs->n_series, cs->capacity, cs->s, cs->e, cs->next);
    fflush(stdout);
    */

    p = cs;
    cs = cs->next;
    if(p->n_series <= 1) {
      p->next = NULL;
      free_clusters(p);
    } else {
      p->next = r_cs;
      r_cs = p;
    }
  }

  free(buf);

  return r_cs;
}

int number_of_clusters(cluster* cs) {
  int n=0;
  for(n=0; cs!=NULL; n++, cs=cs->next)
    ;
  return n;
}

static void most_probable_state_seq(int L, int nh_states, double pbh[],
				    int out_states[]) {
  /* pbh has length nh_states*L, every nh_states is one distribution.
     find the most probable states for the L points and output to out_states.
   */
  int i=0;
  for(i=0; i<L; i++) {
    out_states[i] = arg_max_of(nh_states, pbh + i*nh_states);
  }
}

static void EM_init_distributions(int ns, int n_states[], int nh_states,
				  int np_states,
				  double ph[], double* px[]) {
  /* initialize the distribution randomly */
  int i=0, j=0, k=0;
  double* p=NULL;

  if(ph != NULL) {
    for(i=0; i<np_states; i++) {
      p = ph + i*nh_states;
      for(j=0; j<nh_states; j++) {
	p[j] = 0.25 + 0.5*genrand_real1();
      }
      normalize_dist(nh_states, p);
    }
  }

  if(px != NULL) {
    for(k=0; k<ns; k++) {
      p = px[k];
      for(j=0; j<nh_states; j++) {
	for(i=0; i<n_states[k]; i++) {
	  p[i] = 0.25 + 0.5*genrand_real1();
	}
	normalize_dist(n_states[k], p);
	p += n_states[k];
      }
    }
  }
}

static double EM_E_step(int ns, int n_states[], int nh_states,
			int L, int cs[], int ps[],
			double pbh[], double ph[], double* px[]) {
  /* E-step, calculate pbh, assuming ph and px */
  int i=0, j=0, k=0, x=0;
  double sum=0, tmp=0, log_likelihood = 0;

  for(i=0; i<L; i++) {
    sum = 1e-20; /* to avoid taking log of 0 */
    for(j=0; j<nh_states; j++) {
      tmp = (ps == NULL) ? ph[j] : ph[j + ps[i]*nh_states];
      for(k=0; k<ns; k++) {
	x = cs[i + k*L];
	if(x >= 0) {
	  tmp *= px[k][x + j*n_states[k]];
	}
      }
      pbh[j + i*nh_states] = tmp;
      sum += tmp;
    }
    /* also calculate the log-likelihood */
    log_likelihood += log(sum);

    normalize_dist(nh_states, pbh + i*nh_states);
  }

  return log_likelihood;
}

static void EM_M_step(int ns, int n_states[], int nh_states,
		      int L, int cs[], int np_states, int ps[],
		      double pbh[], double ph[], double* px[]) {
  /* M-step, update distribution. like maximum likelihood, but use pbh
     instead of counts */
  int i=0, j=0, k=0, x=0;
  int* p=NULL;

  if(ps == NULL) np_states=1;

  /* update px */
  for(k=0; k<ns; k++) { /* for each variable */
    memset(px[k], 0, sizeof(double)*n_states[k]*nh_states);
    p = cs + k*L;
    for(i=0; i<L; i++) {
      if(p[i] >= 0) {
	x = p[i];
	for(j=0; j<nh_states; j++) {
	  px[k][x + j*n_states[k]] += pbh[j + i*nh_states];
	}
      }
    }

    for(j=0; j<nh_states; j++) {
      normalize_dist(n_states[k], px[k] + j*n_states[k]);
    }
  }

  /* update ph */
  memset(ph, 0, sizeof(double)*nh_states*np_states);
  if(ps == NULL) {
    for(i=0; i<L; i++) {
      for(j=0; j<nh_states; j++)
	ph[j] += pbh[j + i*nh_states];
    }
  } else {
    for(i=0; i<L; i++) {
      for(j=0; j<nh_states; j++)
	ph[j + ps[i]*nh_states] += pbh[j + i*nh_states];
    }
  }

  for(i=0; i<np_states; i++) normalize_dist(nh_states, ph+i*nh_states);
}

static double perform_EM(int ns, int n_states[], int nh_states,
			 int L, int cs[],
			 int np_states, int ps[],
			 double pbh[], double ph[], double* px[],
			 int n_EM_iter, int n_no_change_restart,
			 int cur_probable_states[], int prev_probable_states[],
			 int out_probable_states[]) {
  /* If h has parent, np_states is the number of states of the
     (grouped) parent, and ps is the values of the (group) parent. And
     ph should have length np_states*nh_states, and will record the
     estimated conditional distribution of h conditional on its
     (grouped) parent.

     If h has no parent, np_states is ignored and ps should be
     NULL. And ph has length nh_states and will record the estimated
     marginal distribution of h.

     refer to re_estimate_val() for the other parameters and their lengths. */
  int iter=0, no_change=0;
  double log_likelihood=0, best_log_likelihood=0;

  /* debug */
  /*
  printf("n_states:"); pr_ints(ns, n_states);
  */

  /* debug */
  /*
  if(ps != NULL) {printf("ps:"); pr_ints(L, ps);}
  */

  for(iter=1; iter<=n_EM_iter; iter++) {
    /* debug */
    /*
    {
      int k=0, j=0;
      printf("** ph: "); pr_doubles(nh_states*np_states, ph);
      for(k=0; k<ns; k++) {
	printf("** px[%d]:\n", k);
	for(j=0; j<nh_states; j++) {
	  printf("    h=%d: ", j);
	  pr_doubles(n_states[k], px[k] + j*n_states[k]);
	}
      }
    }
    */

    log_likelihood = EM_E_step(ns, n_states, nh_states, L, cs, ps, pbh, ph, px);
    printf("EM Iteration %d, log-likelihood: %g\n", iter, log_likelihood);

    /* get the most probable states for h */
    most_probable_state_seq(L, nh_states, pbh, cur_probable_states);

    /* debug */
    /*
    printf("cur_probable_states:"); pr_ints(L, cur_probable_states);
    */

    /* update the best so far, if applicable */
    if(iter==1 || log_likelihood > best_log_likelihood) {
      memcpy(out_probable_states, cur_probable_states, sizeof(int)*L);
      best_log_likelihood = log_likelihood;
    }

    /* check to see if need re-start */
    if(memcmp(cur_probable_states, prev_probable_states, sizeof(int)*L) == 0) {
      no_change++;
    } else {
      no_change = 0;
      memcpy(prev_probable_states, cur_probable_states, sizeof(int)*L);
    }

    if(no_change >= n_no_change_restart) {
      printf("*** most probable states not changed in %d iterations, re-start\n", no_change);
      EM_init_distributions(ns, n_states, nh_states, np_states, ph, px);
    } else {
      /* update ph and px */
      EM_M_step(ns, n_states, nh_states, L, cs, np_states, ps, pbh, ph, px);
    }
  }

  return best_log_likelihood;
}

static void grouped_parents(array2d* data, int n_segs, int lens[],
			    int n_parents, int idx[], int offsets[],
			    int minO, int maxO,
			    int L, int out_states[],
			    int ibuf[]) {
  /* data (in column major form) has expression of n_segs segments,
     with lengths in lens.

     There are n_parents, with indices (with respect to data) in idx,
     and the offsets relative to minO in offsets.

     minO and maxO are the minimum and maximum offsets of the
     children, with respect to the layout in out_states, refer to
     re_estimate_val for their meaning.

     Get the shifted and grouped states of parents and output to
     out_states.  Use 0 if any of the parents is non-overlapping for a
     time points. There are L time points.

     ibuf is for temporary use, should have length 2*n_parents.
   */
  int i=0, j=0, k=0, s=0;
  int eL=0, seg=0, cL=0, sL=0, min_off=0, max_off=0;
  int *s_buf=NULL, *ns_buf=NULL;

  eL = maxO - minO; /* extra length of each segment */
  memset(out_states, 0, sizeof(int)*L);

  min_off = min_int_of(n_parents, offsets);
  max_off = max_int_of(n_parents, offsets);

  s_buf = ibuf;
  ns_buf = ibuf + n_parents;

  for(j=0; j<n_parents; j++) {
    ns_buf[j] = 1 + max_int_of(data->rows, get_col(data, idx[j]));
  }

  /* debug */
  /*
  printf("grouped_parents(): ns_buf:"); pr_ints(n_parents, ns_buf);
  */

  cL=0;
  k=0;
  for(seg=0; seg<n_segs; seg++) {
    sL = eL + lens[seg];
    for(i=0; i<sL; i++) {
      /* i=0 corresponds to minO */
      s = 0;
      if((i-max_off >= 0) &&
	 (i-min_off < lens[seg])) {
	/* get the states of parents */
	for(j=0; j<n_parents; j++) {
	  s_buf[j] = caref(data, cL + i - offsets[j], idx[j]);
	}
	s = states_to_index(n_parents, s_buf, ns_buf);

	/* debug */
	/*
	printf("grouped_parents(): i: %d, s: %d, s_buf:", i, s);
	pr_ints(n_parents, s_buf);
	*/
      }
      /* ** */
      out_states[k + i] = s;
    }
    k += sL;
    cL += lens[seg];
  }
}

static int find_potential_parents(array2d* data, int n_segs, int lens[],
				  int min_offset, int max_offset,
				  int eL, int h[],
				  int max_parents, double score_threshold,
				  links ls_buf[], links* p_links[],
				  int buf1[], int buf2[]) {
  /* for hidden common cause.
     min_offset and max_offset are relative to h.
     eL is extra length of each segment in h.

     Find up to max_parents parents (with offset), where the G2 score
     if >= score_threshold.

     ls_buf and p_links (both length max_parents*(max_offset-min_offset+1)),
     buf1 and buf2 (both length same as h) are for temporary use.

     return the number of potential parents identified, and their info
     are placed in the first few entries of p_links, which point into
     ls_buf.
   */
  int i=0, k=0, ng=0, off=0, tL=0;
  double score=0, test_value=0;

  ng = data->cols;

  k = 0;
  for(i=0; i<ng; i++) {
    /* assume min_offset >= eL */
    for(off=min_offset; off<=max_offset; off++) {
      /* find shifted series */
      tL = shift_by_offset(h, n_segs, lens, off, eL, off, buf1);
      if(tL > 0) {
	shift_by(get_col(data, i), n_segs, lens, 0, off-eL, buf2);
	/* G2 test */
	test_func_gtest(tL, buf1, buf2, &score, &test_value);
	/* fill in link */
	ls_buf[k].from = i;
	ls_buf[k].to = ng+1;
	ls_buf[k].delay = off; /* in fact store the offset */
	ls_buf[k].score = score;
	ls_buf[k].test_value = test_value;
	p_links[k] = ls_buf+k;

	printf(" ** %d -> h, offset: %d\tscore: %g\ttest_value: %g\n",
	       i, off, score, test_value);

	k++;
      }
    }
  }

  qsort(p_links, k, sizeof(links*), decreasing_link_score);

  /* use those with better scores */
  for(i=0; i<max_parents && i<k; i++) {
    if(p_links[i]->score < score_threshold) break;
  }

  return i; /* the number of potential parents */
}

void re_estimate_val(cluster* c, array2d* data, int n_segs, int lens[],
		     int min_delay, int max_delay, int nh_states,
		     int n_EM_iter, int n_no_change_restart,
		     int h_max_parents, double hp_score_threshold,
		     int* out_eL) {
  /* use SVD to do rank-1 approximation to the overlapping parts of
     the series, to get the coefficients for the series.
     buf and buf2 are for tempoary use, and should have length c->cols.

     nh_states is the number of states for the hidden cause.

     output the extra length of each segment to out_eL.
  */
  int i=0, j=0, k=0, ns=0, m=0, np=0;
  int eL=0;
  int minO=0, maxO=0, L=0, seg=0;
  int np_states=0;
  series* se=NULL;
  int *p=NULL, *q=NULL;
  int *ibuf=NULL, *cs=NULL, *n_states=NULL; /* cs and n_states point into ibuf */
  int *cur_probable_states=NULL, *prev_probable_states=NULL; /* point into ibuf */
  int *h_parents_states=NULL; /* point into ibuf */
  int *idx=NULL, *offsets=NULL, *tmp_ibuf=NULL; /* point into ibuf */
  double* buf=NULL;
  double *ph=NULL, *pbh=NULL; /* point into buf, need not release separately */
  double** px=NULL;
  double *cph=NULL;
  double best_log_likelihood=0;

  links *ls_buf=NULL, **p_links=NULL;

  m = data->rows; /* should be equal to sum of lens */
  ns = c->n_series;

  /* assume c->s is the maximum delay allowed. And the original data
     in c is ignored. Now we re-fill it with estimated hidden cause by
     EM, and each segment has length (maxO - minO) + lens[i].

     Assume c has enough capacity, which is true if it is created by
     our other routines.
  */

  printf("Re-estimate values for cluster %d\n", c->id);

  /* minO is min offset, maxO is max offset */
  minO = 2*max_delay;
  maxO = 0;

  for(i=0, se=c->merged; se!=NULL; se=se->next, i++) {
    printf("series %d at offset %d\n", se->idx, se->offset);
    if(se->offset > maxO) maxO = se->offset;
    if(se->offset < minO) minO = se->offset;
  }

  /* also reset c->s and c->e based on the list of merged series */
  c->s = minO; c->e = maxO + m;

  eL = maxO - minO; /* the extra length for each segment */
  if(out_eL) *out_eL = eL;
  L = m + n_segs*eL;

  /* debug */
  /*
  printf("minO: %d, maxO: %d, eL: %d, L: %d\n", minO, maxO, eL, L);
  fflush(stdout);
  */

  /* copy the data */
  j = ns + L*ns + 3*L + 4*h_max_parents;
  ibuf = Malloc(j, int, "re_estimate_val()");
  memset(ibuf, 0, sizeof(int)*j);
  n_states = ibuf; j = ns;
  cs = ibuf + j; j += L*ns;
  cur_probable_states = ibuf + j; j += L;
  prev_probable_states = ibuf + j; j += L;
  h_parents_states = ibuf + j; j += L;
  idx = ibuf + j; j += h_max_parents;
  offsets = ibuf + j; j += h_max_parents;
  tmp_ibuf = ibuf + j; j += 2*h_max_parents;
  /* cs[j + i*L] is for the i-th series at compacted position j.
     -1 means no data, because that part of the series is not overlapping
   */
  p = cs;
  for(k=0, se=c->merged; se!=NULL; se=se->next, k++) {
    q = get_col(data, se->idx);
    /* number of states for the children */
    n_states[k] = 1 + max_int_of(m, q);
    /* segment by segment */
    for(seg=0; seg<n_segs; seg++) {
      for(i=0; i<(se->offset - minO); i++) p[i] = -1;
      for(j=0; j<lens[seg]; i++, j++)
	p[i] = q[j];
      for(; i<(eL + lens[seg]); i++) p[i] = -1;
      p += eL + lens[seg];
      q += lens[seg];
    }
  }

  /* let the hidden cause by h, the children are x_1, x_2, ..., x_ns.
     assume each x_i depends on and only on h, i.e. naive bayes assumption.
   */
  if(nh_states < 2) nh_states = max_int_of(ns, n_states);
  printf("== nh_states: %d\n", nh_states);
  printf("== n_states: "); pr_ints(ns, n_states);

  /* nh_states is number of states for the hidden cause */
  /* ph: marginal of h, length nh_states*np_states
     px[i]: CPT for P(x_i|h)
     pbh[j+t*nh_states]: P(h(t)=j|x_1(t),...,x_ns(t)), for time point t. length L*nh_states.
   */
  /* allocate buffers for distributions */
  px = Malloc(ns, double*, "re_estimate_val()");

  j = nh_states + L*nh_states;
  for(k=0; k<ns; k++) {
    j += n_states[k]*nh_states;
  }

  buf = Malloc(j, double, "re_estimate_val()");
  memset(buf, 0, sizeof(double)*j);

  ph = buf; j = nh_states;
  pbh = buf + j; j += L*nh_states;
  for(k=0; k<ns; k++) {
    px[k] = buf + j;
    j += n_states[k]*nh_states;
  }

  memset(c->val, 0, sizeof(int)*(c->capacity));
  /* use EM to estimate the states of h on the time points, basically
     like naive bayes clustering.
     First without assuming h having any parents.
  */
  EM_init_distributions(ns, n_states, nh_states, 1, ph, px);
  memset(cur_probable_states, 0, sizeof(int)*L);
  memset(prev_probable_states, 0, sizeof(int)*L);

  best_log_likelihood = perform_EM(ns, n_states, nh_states,
				   L, cs, 1, NULL,
				   pbh, ph, px,
				   n_EM_iter, n_no_change_restart,
				   cur_probable_states, prev_probable_states,
				   c->val);
  printf("Best log-likelihood found (no parents): %g\n", best_log_likelihood);
  /* find potential parents for h
   */
  printf("== Try to find potential parents for the estimated hidden cause.\n");

  j = (data->cols)*(max_delay + 1);
  ls_buf = Malloc(j, links, "re_estimate_val() ls_buf");
  p_links = Malloc(j, links*, "re_estimate_val() p_links");

  k = eL < max_delay ? max_delay : eL;
  np = find_potential_parents(data, n_segs, lens,
			      k + min_delay, k + max_delay,
			      eL, c->val, h_max_parents, hp_score_threshold,
			      ls_buf, p_links,
			      cur_probable_states, prev_probable_states);

  printf("== Found %d potential parents.\n", np);
  /* if found, use EM again */
  if(np > 0) {
    /* extract the indices and offsets */
    for(i=0; i<np; i++) {
      idx[i] = p_links[i]->from;
      offsets[i] = p_links[i]->delay;
    }
    printf("They are:"); pr_ints(np, idx);

    printf("Use EM to re-estimate the hidden cause with potential parents.\n");

    /* group the parent states */
    grouped_parents(data, n_segs, lens, np, idx, offsets,
		    minO, maxO, L, h_parents_states,
		    tmp_ibuf);
    /* EM again, re-initialize ph, but keep px */
    np_states = 1+max_int_of(L, h_parents_states);

    /* debug */
    /*
    printf(" ** np_states: %d\tnh_states: %d\tL: %d\n", np_states, nh_states, L);
    */

    cph = Malloc(np_states*nh_states, double, "re_estimate_val() cph");

    EM_init_distributions(ns, n_states, nh_states, np_states, cph, px);
    memset(cur_probable_states, 0, sizeof(int)*L);
    memset(prev_probable_states, 0, sizeof(int)*L);

    best_log_likelihood = perform_EM(ns, n_states, nh_states,
				     L, cs, np_states, h_parents_states,
				     pbh, cph, px,
				     n_EM_iter, n_no_change_restart,
				     cur_probable_states, prev_probable_states,
				     c->val);
    printf("Best log-likelihood found (with parents): %g\n",
	   best_log_likelihood);
    free(cph);
  }
  /* cleanup */
  free(ibuf);
  free(buf);
  free(px);
  free(ls_buf);
  free(p_links);
}

array2d* augment_data_with_clusters(array2d* data, int n_segs, int lens[],
				    int candidates[], int n_candidates, int n_candidates_parents,
				    int min_delay, int max_delay, int nh_states,
				    int n_EM_iter, int n_no_change_restart,
				    double threshold,
				    int h_max_parents, double hp_score_threshold,
				    char** out_forbidden) {
  /* data is assumed to be in column major form, with n_segs segments,
     and the lengths of the segments are in lens.

     Consider the first n_candidates column indices in candidates for clustering.
     return a possibly augmented data (still in column major form),
     that adds the cluster centers as columns.

     May return the original data if no columns are added.

     nh_states is the number of states for the hidden cause.

     threshold is for testing if the association score is large enough to include to a cluster.

     If out_forbidden is non-NULL, *out_forbidden will be written a newly allocated array
       where out_forbidden[i*ng + j] is 1 if i->j is forbidden,
       where ng is the number of columns in augmented data.
       The array in *out_forbidden should be freed after use.
  */
  int n_cs=0, m=0, n=0, i=0, s=0;
  array2d* a_data=NULL;
  cluster *cs=NULL, *p=NULL;
  series* s_buf=NULL;
  int seg=0, cL=0;
  int *u=NULL, *v=NULL, extra_L=0;

  m = data->rows;
  n = n_candidates + n_candidates_parents;

  s_buf = Malloc(n, series, "augment_data_with_clusters()");

  cs = find_clusters(data, n_segs,lens, candidates, n_candidates,
		     max_delay, threshold, s_buf);
  n_cs = number_of_clusters(cs);

  if(n_cs > 0) {
    a_data = new_array2d(m, data->cols + n_cs); /* n_cs new columns */

    /* copy the original columns */
    for(i=0; i<(data->cols); i++) {
      memcpy(get_col(a_data,i), get_col(data,i), sizeof(int)*m);
    }

    /* the new columns */
    for(p=cs, i=(data->cols); p!=NULL; p=p->next, i++) {
      printf("Add cluster %d as column %d\n", p->id, i);
      memset(get_col(a_data, i), 0, sizeof(int)*m);

      /* re-estimate the cluster values from the series */
      re_estimate_val(p, data, n_segs,lens, min_delay, max_delay, nh_states,
		      n_EM_iter, n_no_change_restart,
		      h_max_parents, hp_score_threshold,
		      &extra_L);

      /* for the time being, simply use the last part of the merged
	 series for each segment, so that it is before the merged
	 series, but with at least one more delay.
      */
      u = get_col(a_data, i);
      v = p->val;
      for(seg=0, cL=0; seg<n_segs; seg++) {
	s = extra_L + 1; /* shift at least one */
	if(s < max_delay) s = max_delay;
	printf("Taking from %d to %d, length %d.\n", cL+s, cL+s+lens[seg]-1, lens[seg]-1);
	memcpy(u, v + s, sizeof(int)*(lens[seg]-1));

	/* next segment */
	cL += lens[seg] + extra_L;
	u += lens[seg];
	v += lens[seg] + extra_L;
      }
    }
  } else {
    a_data = data;
  }

  /* record the forbidden links, because they are likely children of
     an introduced hidden common node */
  if(out_forbidden) {
    {
      char* f=NULL;
      int ng = a_data->cols;
      series *a=NULL, *b=NULL;
      f = Malloc(ng*ng, char, "augment_data_with_clusters() out_forbidden");
      memset(f, 0, ng*ng*sizeof(char));
      for(p=cs, i=(data->cols); p!=NULL; p=p->next, i++) {
	/* brute force */
	for(a=p->merged; a!=NULL; a=a->next) {
	  for(b=a->next; b!=NULL; b=b->next) {
	    f[ng*(a->idx) + (b->idx)] = 1; /* a -> b forbidden */
	    f[ng*(b->idx) + (a->idx)] = 1; /* b -> a forbidden */
	    printf("*** forbidden links %d <--> %d\n", a->idx, b->idx);
	  }
	  f[ng*(a->idx) + i] = 1; /* children back to the hidden common cause forbidden */
	}
      }
      *out_forbidden = f;
    }
  }

  /* done, cleanup */
  free_clusters(cs);
  free(s_buf);

  return a_data;
}

int hcc_dclinde(array2d* d, int n_segs, int lens[],
		int min_delay, int max_delay, int nh_states,
		int n_EM_iter, int n_no_change_restart,
		int method1, int method2,
		double st1, double st2, int max_n, int pruning,
		int is_one_delay, int is_no_dup, int is_keep_zero_delay,
		int i_start, int i_end, double alpha, int allow_self_link,
		int is_init_only, int is_no_cpt,
		double eV, double eV_tolerance, double eV_threshold,
		int h_max_parents, double hp_score_threshold,
		links** out_links,
		char* out_grn, char* out_init_grn) {
  /* main function of HCC-DCLINDE.

     d is n by g array, the n_segs segments of expression data in column major format,
     where rows are time points, and columns are genes.
     The first lens[0] rows are segment 1,
     the next lens[1] rows are segment 2, and so on.

     If eV > 0, it is the expected bias of the conditional
     distributions, and if the bias of a gene is different from
     expected, i.e. |eV - bias| > eV_tolerance, it will be candidate
     of having hidden common node.

     If eV <= 0, the expected bias of the error terms will be
     estimated as the median.

     eV_threshold is for testing if the association score is large
     enough to include to a cluster.

     nh_states is the number of states for the hidden cause.

     Other parameters are similar to dclinde() above.

     Some parameters are for D-CLINDE, some are for GlobalMIT+ / GlobalMIT*.

     To try to handle hidden common nodes in GRN.
     Main idea:
     First use DCLINDE to infer an GRN, then check the error levels of each gene.

     Those genes with error levels very different from expected, and
     their parents which have no parents are considered candidates for
     having hidden common nodes and the links with them are cut.

     Then cluster the series of the candidates to estimate the hidden
     common node(s). And then run CLINDE on these candidates and
     estimated hidden common nodes, with some prunings.

     Finally merge the two GRNs to give the final GRN.

     If out_links is non-NULL, allocate a copy of the final links to
     *out_links.
   */
  int ng=0;
  int n_links=0, i=0, j=0, k=0;
  links* the_links=NULL;

  double *bias_var=NULL;
  node* grn_ns=NULL;
  int n_candidates=0, n_candidates_parents=0, *candidates=NULL, *is_in=NULL;
  int ns_links=0; /* for the sub part involving hidden common node */
  links* s_links=NULL;
  int n_final_links=0;
  array2d* a_data=NULL;

  ng = d->cols; /* number of observed genes */
  a_data = d; /* may later change to point to augmented data */

  /* get the inititial network */
  printf("==== Initial GRN\n");
  if(method1 == METHOD_DCLINDE) {
    n_links = dclinde(d, n_segs, lens, st1, st2, min_delay, max_delay, max_n,
		      pruning, is_one_delay, is_no_dup, is_keep_zero_delay,
		      NULL, &the_links, is_init_only ? (out_init_grn==NULL ? out_grn : out_init_grn) : out_init_grn, NULL);
  } else if(method1 == METHOD_GLOBALMIT_FAST || METHOD_GLOBALMIT) {
    n_links = globalMIT_exe(d, n_segs, lens, max_delay, i_start, i_end, alpha, allow_self_link,
			    method1 == METHOD_GLOBALMIT_FAST ? globalMIT_fixedOrder_exe_wrapper : globalMIT_exe_wrapper,
			    NULL, &the_links,
			    is_init_only ? (out_init_grn==NULL ? out_grn : out_init_grn) : out_init_grn);
  }
  printf("==== Learnt Initial GRN\n");
  print_n_links(stdout, n_links, the_links);
  printf("==== End Initial GRN\n");

  /* the CPT, also grn_ns records the links in a different way */
  grn_ns = links_to_nodes(ng, n_links, the_links, d, n_segs, lens);
  if(!is_no_cpt) {
    printf("==== CPT for Initial GRN\n");
    print_nodes(stdout, ng, grn_ns);
    printf("==== End CPT for Initial GRN\n");
  }

  /* estimate the bias */
  {
    int mL=1;
    double* tmp_buf=NULL;
    printf("=== Estimated Bias\n");
    bias_var = Malloc(2*ng, double, "hcc_dclinde()");

    for(i=0; i<ng; i++) {
      if(grn_ns[i].n_cpt > mL) mL = grn_ns[i].n_cpt;
    }

    tmp_buf = Malloc(mL, double, "hcc_dclinde() estimate bias");

    for(i=0; i<ng; i++) {
      bias_var[i] = estimate_bias(grn_ns[i].n_states,
				  grn_ns[i].n_cpt,
				  grn_ns[i].cpt,
				  tmp_buf);
      bias_var[i+ng] = bias_var[i];
      printf("Gene %d: %g\n", i, bias_var[i]);
    }
    printf("=== End Estimated Bias\n");

    free(tmp_buf);
  }

  /* threshold for bias */
  if(eV <= 0) { /* use the median as estimate */
    qsort(bias_var+ng, ng, sizeof(double), decreasing_double);
    eV = bias_var[ng+(ng/2)];
  }
  printf("Bias thresholds: %g +- %g = [%g , %g]\n",
	 eV, eV_tolerance, eV - eV_tolerance, eV + eV_tolerance);

  /* determine candidate genes as having hidden common node.
     First mark them, then record their index.
   */
  n_candidates = 0;
  candidates = Malloc(2*ng, int, "hcc_dclinde()");
  memset(candidates, 0, 2*ng*sizeof(int));
  is_in = candidates + ng; /* is_in need not be freed separately */
  for(i=0; i<ng; i++) {
    if((bias_var[i] < eV-eV_tolerance) ||
       (bias_var[i] > eV+eV_tolerance)) {
      is_in[i] = 1;
      n_candidates++;
    }
  }
  n_candidates_parents = 0;
  /* also the parents of the candidates (which are not candidates) are
     also included, but restricted */
  for(i=0; i<ng; i++) {
    if(is_in[i] == 1) {
      for(j=grn_ns[i].n_parents-1;  j>=0; j--) {
	k = grn_ns[i].parents[j];
	/* k --> i */
	if(is_in[k]==0) {
	  is_in[k] = 2;
	  n_candidates_parents++;
	}
      }
    }
  }

  j=0;
  for(i=0; i<ng; i++) {if(is_in[i]==1) candidates[j++] = i;}
  for(i=0; i<ng; i++) {if(is_in[i]==2) candidates[j++] = i;}

  printf("==== Unexpected bias in %d gene(s):", n_candidates);
  for(i=0; i<n_candidates; i++) printf(" %d", candidates[i]);
  printf("\n");
  printf("==== Their parents, %d gene(s):", n_candidates_parents);
  for(i=0; i<n_candidates_parents; i++)
    printf(" %d", candidates[i+n_candidates]);
  printf("\n");

  report_stop_watch();

  /* do the hidden common node part, if applicable */
  n_final_links = n_links;
  if(n_candidates > 1 && !is_init_only) {
    {
      char* forbidden=NULL;
      int n=0, n_ng=0;

      n = n_candidates + n_candidates_parents;
      printf("====== Clustering and Augmenting the Data ======\n");
      a_data = augment_data_with_clusters(d, n_segs,lens,
					  candidates, n_candidates,
					  n_candidates_parents,
					  min_delay, max_delay, nh_states,
					  n_EM_iter, n_no_change_restart,
					  eV_threshold,
					  h_max_parents, hp_score_threshold,
					  &forbidden);

      report_stop_watch();
      n_ng = a_data->cols;
      if(n_ng > n) { /* introduced hidden common node, re-learn the network */
	printf("==== GRN for potential hidden common causes\n");
	if(method2 == METHOD_DCLINDE) {
	  ns_links = dclinde(a_data, n_segs, lens, st1, st2, min_delay, max_delay, max_n,
			     pruning, is_one_delay, is_no_dup, is_keep_zero_delay,
			     forbidden, &s_links, NULL, NULL);
	} else if(method2 == METHOD_GLOBALMIT_FAST || METHOD_GLOBALMIT) {
	  ns_links = globalMIT_exe(a_data, n_segs, lens, max_delay, 0, n_ng-1, alpha, allow_self_link,
				   method2 == METHOD_GLOBALMIT_FAST ? globalMIT_fixedOrder_exe_wrapper : globalMIT_exe_wrapper,
				   forbidden, &s_links, NULL);
	}
	printf("==== End GRN for potential hidden common causes\n");
	report_stop_watch();

	/* replace the original links */
	free(the_links);
	n_final_links = ns_links;
	the_links = s_links; /* transfer ownership */
	s_links = NULL;
      } else {
	printf("==== No candidate hidden common node identified.\n");
	free(s_links);
      }
      /* small clean up */
      free(forbidden);
    }
  }

  /* output */
  printf("==== Final GRN\n");
  for(i=0; i<n_final_links; i++) {
    print_link(stdout, the_links+i);
  }
  printf("==== End Final GRN\n");

  if(out_grn != NULL) {
    {
      FILE* tmpf=NULL;
      printf("=== output final GRN to file %s.\n", out_grn);
      tmpf = fopen(out_grn, "wt");
      if(tmpf == NULL) {
	fprintf(stderr, "Error: Cannot open file %s to output GRN.\n",
		out_grn);
      } else {
	for(i=0; i<n_final_links; i++) {
	  print_link(tmpf, the_links+i);
	}
	fclose(tmpf);
      }
    }
  }

  if(!is_no_cpt) {
    printf("==== CPT for Final GRN\n");
    calc_and_print_CPT(a_data->cols, n_final_links, the_links,
		       a_data, n_segs, lens);
    printf("==== End CPT for Final GRN\n");
  }

  if(out_links!=NULL && n_final_links>0) {
    {
      links* outL = NULL;
      int k=0;
      outL = Malloc(n_final_links, links, "hcc_clinde()");

      for(i=0; i<n_final_links; i++) {
	outL[k++] = the_links[i];
      }
      *out_links = outL;
    }
  }

  /* clean up */
  free(the_links);
  free(bias_var);
  free(candidates);

  if(a_data != d) free_array2d(a_data);
  free_nodes(ng, grn_ns);

  /* done */
  return n_final_links;
}

/***************************************************************/
/* discrete data in form of small integers */

array2d* read_expression_data(char* filename) {
  /* return an 2d array (row major) if can read properly. return NULL on error. */
  array2d* r=NULL;
  int rows=0, cols=0;
  ELEM* res=read_dTSV(filename, &rows, &cols);
  if(res == NULL) {
    fprintf(stderr, "Error: Cannot properly read %s for expression data.\n", filename);
    return NULL;
  }
  r = Malloc(1, array2d, "read_expression_data()");
  r->rows = rows;
  r->cols = cols;
  r->v = res;
  return r;
}

array2d* read_n_data(int n, char* fn[], int out_lens[]) {
  /* read n expression data, fn are the filenames.
     The expression data should have the same number of columns.
     Put them one by one in the same 2D array, and output the number
     of rows of each data in out_lens.
     Return the array in column major format, so that it is easier to extract a column.
     Return the 2D array if OK. Return NULL if error.
   */
  int i=0,j=0,k=0, ncols=0, nrows=0, nr=0;
  array2d* r = NULL;
  array2d** buf = NULL;
  buf = Malloc(n,array2d*, "read_n_data()");
  memset(buf, 0, sizeof(array2d*)*n);

  printf("About to read %d expression data\n", n);

  for(i=0; i<n; i++) {
    printf("Reading the %dth expression data %s ... ", i+1, fn[i]);
    buf[i] = read_expression_data(fn[i]);
    if(buf[i] == NULL) goto done;
    if(i==0) {
      ncols = buf[i]->cols;
    } else if(ncols != buf[i]->cols) {
      fprintf(stderr, "Error: the %dth expression data %s has %d columns, but previous ones have %d columns.\n",
	      (i+1),fn[i], buf[i]->cols, ncols);
      goto done;
    }
    printf("OK.\n  Read %d rows, %d cols.\n", buf[i]->rows, buf[i]->cols);
    nrows += buf[i]->rows;
  }
  /* copy them to one array, in column major format */
  r = new_array2d(nrows, ncols);
  for(i=0, nr=0; i<n; i++) {
    out_lens[i] = buf[i]->rows;;
    for(j=0; j<buf[i]->rows; j++, nr++)
      for(k=0; k<ncols; k++)
	caset(r, nr,k, aref(buf[i], j,k));
  }
  printf("Done\nRead %d time points, %d genes.\n", nrows, ncols);
  /***/
 done:
  for(i=0; i<n; i++)
    if(buf[i]) free_array2d(buf[i]);
  free(buf);
  return r;
}

/***************************************************************/
int increasing_int(const void* a,const void* b) {
  return *((int*)a) - *((int*)b);
}

static int remove_neighbor_dups(int L, int x[]) {
  /* if two neighboring element in x are the same, retain only one.
     compact x in place.
     return the new length.
   */
  int i=0, j=0;

  if(L <= 0) return 0;

  for(i=1; i<L; i++) {
    if(x[i] != x[j]) {
      x[++j] = x[i];
    }
  }

  return j+1;
}

static int position_of(int a, int L, int x[]) {
  int i=0;
  for(i=0; i<L; i++)
    if(x[i] == a) return i;
  return L;
}

static int rename_states(int L, int x[]) {
  /* replace the states in x to be small non-negative integers, if applicable.
     returns the number of states.
   */
  int* p=NULL;
  int i=0, n=0;

  if(L <= 0) return 0;

  p = reserve_temp_buffer(L);
  memcpy(p, x, sizeof(int)*L);
  qsort(p, L, sizeof(int), increasing_int);

  n = remove_neighbor_dups(L, p);
  if(p[0]==0 && p[n-1]==n-1) return n; /* no need to rename */

  printf("*** Rename");
  for(i=0; i<n; i++) printf(", %d=>%d", p[i],i);
  printf(" ***");
  for(i=0; i<L; i++)
    x[i] = position_of(x[i], n,p);

  return n;
}

static void rename_data_states(array2d* d) {
  /* assume d is in column major format */
  int i=0, n=0, s=0;
  n = d->cols;
  for(i=0; i<n; i++) {
    printf("Gene %d ", i);
    s = rename_states(d->rows, get_col(d,i));
    printf(": %d states\n", s);
  }
}

/***************************************************************/
#define DEFAULT_SEED       "0"
#define DEFAULT_ST1        2
#define DEFAULT_ST2        2
#define DEFAULT_MIN_DELAY  1
#define DEFAULT_MAX_DELAY  4
#define DEFAULT_MAX_N      2
#define DEFAULT_PRUNING    "all"
#define DEFAULT_METHOD1    "mit_fast"
#define DEFAULT_METHOD2    "mit_fast"
#define DEFAULT_ALPHA      0.999
#define DEFAULT_EV         0.65
#define DEFAULT_EV_TOLERANCE 0.05
#define DEFAULT_EV_THRESHOLD 2.3
#define DEFAULT_NH_STATES  0
#define DEFAULT_N_EM_ITER  100
#define DEFAULT_N_NO_CHANGE_RESTART 3
#define DEFAULT_MP         3
#define DEFAULT_HP_ST      2

#define xstr(s) str(s)
#define str(s) #s

/* The option structure of this program */
option all_options[] = {
  /* str,type(STRING),is_required(0),description,index(0),end_index(0), val_int(0),val_real(0),val_str(NULL) */
  {"-?",         FLAG,   0,  "Showing the usage.\n", 0,0,  0,  0.0,    NULL},

  {"-data",      LIST,   1,  "File name(s) of the input expression data (tab/space separated) small non-negative integers for discrete data, each row is a time point, each column is a gene. Each file should have the same number of columns, but may have different number of rows.\n", 0,-1,  -1,  0.0,    NULL},

  {"-seed",      STRING, 0,  "The seed of the pseudo random number generator for EM initial distribution. If contains '-' (.e.g negative number), use current time. Default " xstr(DEFAULT_SEED) "\n", 0,0,  0,  0.0,    DEFAULT_SEED},

  {"-st1",       REAL,   0,  "Score threshold for stage 1 of D-CLINDE. Score is -log10(p), where p is the intended p-value. Default " xstr(DEFAULT_ST1) ".\n", 0,0,  0,  DEFAULT_ST1,    NULL},
  {"-st2",       REAL,   0,  "Score threshold for stage 2 of D-CLINDE. Similar to st1. Default " xstr(DEFAULT_ST2) ".\n", 0,0,  0,  DEFAULT_ST2,    NULL},
  {"-min.delay", INT,    0,  "Minimum delay (lags) to use in the inference (not applicable to GlobalMIT+/*), default " xstr(DEFAULT_MIN_DELAY) ".\n", 0,0,  DEFAULT_MIN_DELAY,  0.0,    NULL},

  {"-max.delay", INT,    0,  "Maximum delay (lags) to use in the inference, default " xstr(DEFAULT_MAX_DELAY) ".\n", 0,0,  DEFAULT_MAX_DELAY,  0.0,    NULL},

  {"-max.n",     INT,    0,  "Maximum number of parents to condition on in stage 2 of D-CLINDE, default " xstr(DEFAULT_MAX_N) ".\n", 0,0,  DEFAULT_MAX_N,  0.0,    NULL},

  {"-pruning",   STRING, 0,  "Pruning strategy in stage 2 of D-CLINDE, can be all (consider all neighbors of the two vertices when doing conditional test of the link between the two) or common (consider only common neighbors of the two vertices). Default " xstr(DEFAULT_PRUNING) "\n", 0,0,  0,  0.0,    DEFAULT_PRUNING},
  /*
  {"-local.delay", FLAG,   0,  "In stage 1, in addition to score threshold, also use local maximum of absolute correlation to filter possible delays. Default false.\n", 0,0,  0,  0.0,    NULL},
  */
  {"-keep.zero.delay", FLAG,  0,  "To keep links with zero delay after stage 2 of D-CLINDE and before remove dups. Default false.\n", 0,0,  0,  0.0,    NULL},
  {"-one.delay", FLAG,   0,  "To keep only one delay (the one with best score, smallest delay) for each link after stage 1 of D-CLINDE. Default false.\n", 0,0,  0,  0.0,    NULL},
  {"-no.dup",    FLAG,   0,  "Remove duplicate links. To keep only one delay (the one with best score, smallest delay) for each link after stage 2 of D-CLINDE. Default false.\n", 0,0,  0,  0.0,    NULL},

  {"-istart",    INT,    0,  "The starting 0-based index of gene for which to find parents. Only for GlobalMIT+/* in initial GRN, default 0.\n", 0,0,  0,  0.0,    NULL},
  {"-iend",      INT,    0,  "The ending (inclusive) 0-based index of gene for which to find parents. Only for GlobalMIT+/* in initial GRN, default -1, which means n_genes-1.\n", 0,0,  -1,  0.0,    NULL},
  {"-alpha",     REAL,   0,  "Significance level for GlobalMIT+/*. Default " xstr(DEFAULT_ALPHA) ".\n", 0,0,  0,  DEFAULT_ALPHA,    NULL},
  {"-no_self_link", FLAG,  0,  "To disallow self link (with delay) for GlobalMIT+/*, default false.\n", 0,0,  0,  0.0,    NULL},

  {"-no.cpt",    FLAG,   0,  "Do not print the conditional distributions for the final GRN. Default false.\n", 0,0,  0,  0.0,    NULL},

  {"-out.grn",   STRING, 0,  "Name of file to output the final GRN. Default NULL, which means do not output to file.\n", 0,0,  0,  0.0,    NULL},
  {"-out.init.grn",STRING, 0,  "Name of file to output the init GRN. Default NULL, which means do not output to file.\n", 0,0,  0,  0.0,    NULL},

  {"-init.only", FLAG,   0,  "To learn only the initial GRN, which does not consider possible hidden common causes. Default false.\n", 0,0,  0,  0.0,    NULL},

  {"-method1",   STRING, 0,  "Method for learning initial GRN, can be either " xstr(M_DCLINDE_NAME) " for D-CLINDE, " xstr(M_MIT_NAME) " for GlobalMIT+, or " xstr(M_MIT_FAST_NAME) " for GlobalMIT*. Default " xstr(DEFAULT_METHOD1) "\n", 0,0,  0,  0.0,    DEFAULT_METHOD1},
  {"-method2",   STRING, 0,  "Method for learning sub-GRN after potential hidden common cause(s) are identified and estimated, can be either " xstr(M_DCLINDE_NAME) " for D-CLINDE, " xstr(M_MIT_NAME) " for GlobalMIT+, or " xstr(M_MIT_FAST_NAME) " for GlobalMIT*. Default " xstr(DEFAULT_METHOD2) "\n", 0,0,  0,  0.0,    DEFAULT_METHOD2},

  {"-eV",        REAL,   0,  "Expected Bias if > 0. If <= 0, the expected bias will be estimated to be the median of the bias of the initial GRN. If bias of a gene is > eV + eV.tolerance or < eV - eV.tolerance, that gene will be considered for having hidden common node. Default " xstr(DEFAULT_EV) ".\n", 0,0,  0,  DEFAULT_EV,    NULL},
  {"-eV.tolerance",        REAL,   0,  "If bias of a gene is > eV + eV.tolerance or < eV - eV.tolerance, that gene will be considered for having hidden common node. Default " xstr(DEFAULT_EV_TOLERANCE) ".\n", 0,0,  0,  DEFAULT_EV_TOLERANCE,    NULL},
  {"-eV.threshold",        REAL,   0,  "Score threshold of association for clustering using G2 test. Association score is -log_10(p-value). Similar to st1 and st2. Default " xstr(DEFAULT_EV_THRESHOLD) ".\n", 0,0,  0,  DEFAULT_EV_THRESHOLD,    NULL},

  {"-nh.states", INT,    0,  "Number of states for hidden cause(s). If < 2, use the maximum number of states of children. Default " xstr(DEFAULT_NH_STATES) ".\n", 0,0,  DEFAULT_NH_STATES,  0.0,    NULL},

  {"-n.EM.iter", INT,    0,  "Number of iterations for EM in estimating hidden cause(s). Default " xstr(DEFAULT_N_EM_ITER) ".\n", 0,0,  DEFAULT_N_EM_ITER,  0.0,    NULL},
  {"-n.no.change.restart", INT,    0,  "In EM, if the most probable hidden states do not change after n.no.change.restart iterations, restart EM. Default " xstr(DEFAULT_N_NO_CHANGE_RESTART) ".\n", 0,0,  DEFAULT_N_NO_CHANGE_RESTART,  0.0,    NULL},

  {"-mp",        INT,    0,  "The maximum number of potential parents to find for hidden cause after the first EM. Default " xstr(DEFAULT_MP) ".\n", 0,0,  DEFAULT_MP,  0.0,    NULL},
  {"-hp.st",     REAL,   0,  "Score threshold of association using G2 test for finding potential parent for hidden cause. Association score is -log_10(p-value). Similar to st1 and st2. Default " xstr(DEFAULT_HP_ST) ".\n", 0,0,  0,  DEFAULT_HP_ST,    NULL},

  {NULL,      STRING, 0,  NULL, 0,0,  0,  0.0,    NULL} /* end marker */
};

/***********/
int our_stricmp(char* a, char* b) {
  int r, i=0;
  do {
    r = tolower(a[i]) - tolower(b[i]);
    if(r > 0) return 1;
    if(r < 0) return -1;
    i++;
  } while(a[i]!='\0' || b[i]!='\0');
  return 0;
}

int match_option(char* str, int n, char* choices[],
		 int warning_if_not_found, FILE* outs, char* option_name) {
  /* return -1 if not found. return the index if found */
  /* option_name is for the warning */
  int i=0;
  for(i=0; i<n; i++) {
    if(our_stricmp(str, choices[i])==0) return i;
  }
  if(!warning_if_not_found) return -1;

  /* not found, indicate the choices */
  fprintf(outs, "Warning: Unknown option \"%s\", possible choices are:\n", str);
  for(i=0; i<n; i++)
    fprintf(outs, "\t%s\n", choices[i]);
  fprintf(outs, "\n");
  return -1;
}

/***************************************************************/
void no_abort_handler (const char * reason, 
		       const char * file, 
		       int line, 
		       int gsl_errno) {
  /* for GSL handler, to print something but do not abort. */
  fprintf(stderr, "GSL: %s:%d: ERROR: (%d) %s\n", file, line, gsl_errno, reason);
  fprintf(stdout, "GSL: %s:%d: ERROR: (%d) %s\n", file, line, gsl_errno, reason);
}

/***************************************************************/

int main(int argc, char* argv[]) {
  double st1=DEFAULT_ST1, st2=DEFAULT_ST2;
  double eV=DEFAULT_EV, eV_tolerance=DEFAULT_EV_TOLERANCE, eV_threshold=DEFAULT_EV_THRESHOLD;
  int h_max_parents = DEFAULT_MP;
  double hp_score_threshold = DEFAULT_HP_ST;
  int nh_states=DEFAULT_NH_STATES;
  int n_EM_iter=DEFAULT_N_EM_ITER, n_no_change_restart=DEFAULT_N_NO_CHANGE_RESTART;
  int min_delay=DEFAULT_MIN_DELAY, max_delay=DEFAULT_MAX_DELAY, max_n=DEFAULT_MAX_N;
  char* s_seed=DEFAULT_SEED;
  char *s_pruning=DEFAULT_PRUNING, *s_method1=DEFAULT_METHOD1, *s_method2=DEFAULT_METHOD2;
  char *s_out_grn=NULL, *s_out_init_grn=NULL;
  int pruning=0, method1=0, method2=0;
  int is_one_delay=0, is_no_dup=0, is_keep_zero_delay=0, is_no_cpt=0, is_init_only=0;
  unsigned long seed=0;

  int istart = 0, iend = -1, is_no_self_link=0;
  double alpha = DEFAULT_ALPHA;

  links* out_links=NULL;

  int s_idx=0, e_idx=-1;
  int n_data=0, *data_lens=NULL;
  array2d *data=NULL;

  /* do some initialization */
  if(parse_options(argc,argv,all_options)) {
    /* error parsing options */
    usage(stderr, argv[0], all_options);
    exit(-1);
  }
  if(get_flag_option_value("-?",all_options,0))
    usage(stdout, argv[0], all_options);

  if(!get_list_option_range("-data", all_options, &s_idx, &e_idx) || (e_idx < s_idx)) {
    fprintf(stderr, "Input expression file names not given.\n");
    usage(stderr, argv[0], all_options);
    exit(-1);
  }
  n_data = e_idx - s_idx; /* both inclusive, but start includes the -data */

  s_seed = get_str_option_value("-seed", all_options, s_seed);

  st1 = get_real_option_value("-st1", all_options, st1);
  st2 = get_real_option_value("-st2", all_options, st2);

  eV = get_real_option_value("-eV", all_options, eV);
  eV_tolerance = get_real_option_value("-eV.tolerance", all_options, eV_tolerance);
  eV_threshold = get_real_option_value("-eV.threshold", all_options, eV_threshold);

  h_max_parents = get_int_option_value("-mp", all_options, h_max_parents);
  hp_score_threshold = get_real_option_value("-hp.st", all_options, hp_score_threshold);

  nh_states = get_int_option_value("-nh.states", all_options, nh_states);
  n_EM_iter = get_int_option_value("-n.EM.iter", all_options, n_EM_iter);
  n_no_change_restart = get_int_option_value("-n.no.change.restart", all_options, n_no_change_restart);

  min_delay = get_int_option_value("-min.delay", all_options, min_delay);
  max_delay = get_int_option_value("-max.delay", all_options, max_delay);
  max_n = get_int_option_value("-max.n", all_options, max_n);
  s_pruning = get_str_option_value("-pruning", all_options, s_pruning);
  s_method1 = get_str_option_value("-method1", all_options, s_method1);
  s_method2 = get_str_option_value("-method2", all_options, s_method2);
  is_one_delay = get_flag_option_value("-one.delay", all_options, is_one_delay);
  is_no_dup = get_flag_option_value("-no.dup", all_options, is_no_dup);
  is_keep_zero_delay = get_flag_option_value("-keep.zero.delay", all_options, is_keep_zero_delay);

  istart = get_int_option_value("-istart", all_options, istart);
  iend = get_int_option_value("-iend", all_options, iend);
  alpha = get_real_option_value("-alpha", all_options, alpha);
  is_no_self_link = get_flag_option_value("-no_self_link", all_options, is_no_self_link);

  is_no_cpt = get_flag_option_value("-no.cpt", all_options, is_no_cpt);
  is_init_only = get_flag_option_value("-init.only", all_options, is_init_only);

  s_out_grn = get_str_option_value("-out.grn", all_options, s_out_grn);
  s_out_init_grn = get_str_option_value("-out.init.grn", all_options, s_out_init_grn);

  if(strchr(s_seed, '-') != NULL) {
    printf("Use current time as seed.\n");
    seed = time(NULL);
  } else {
    seed = atol(s_seed);
  }

  if((pruning = match_option(s_pruning, N_PRUNINGS, pruning_names, 1,stdout,"-pruning")) < 0) {
    printf("Use default for pruning instead.\n");
    pruning = 0;
  }

  if((method1 = match_option(s_method1, N_METHODS, method_names, 1,stdout,"-method1")) < 0) {
    printf("Use default for method1 instead.\n");
    method1 = 0;
  }

  if((method2 = match_option(s_method2, N_METHODS, method_names, 1,stdout,"-method2")) < 0) {
    printf("Use default for method2 instead.\n");
    method2 = 0;
  }

  if(min_delay < 0) {
    printf("min.delay should be non-negative integer, but got %d\n", min_delay);
    printf("Use default instead.\n");
    min_delay = DEFAULT_MIN_DELAY;
  }
  if(max_delay <= 0) {
    printf("max.delay should be positive integer, but got %d\n", max_delay);
    printf("Use default instead.\n");
    max_delay = DEFAULT_MAX_DELAY;
  }
  if(min_delay > max_delay) {
    printf("min.delay should not be larger than max.delay, but got %d and %d\n", min_delay, max_delay);
    printf("Use default for min.delay instead.\n");
    min_delay = DEFAULT_MIN_DELAY;
  }

  /*********/
  start_stop_watch();
  printf("====== HCC-DCLINDE Grn Inference ======\n");
  printf("Seed for pseudo random number generator: %lu\n", seed);
  if(eV > 0)
    printf("Expected Bias: %f\n", eV);
  else
    printf("Expected Bias: Estimated\n");
  printf("Bias Tolerance: %f\n", eV_tolerance);
  printf("Score threshold for clustering: %g\n", eV_threshold);

  printf("Method to learn initial GRN: %s (%s)\n", method_long_names[method1], method_names[method1]);
  if(!is_init_only) {
    printf("Method to learn sub-GRN after identifying and estimating hidden common cause(s): %s (%s)\n",
	   method_long_names[method2], method_names[method2]);
  }

  if(method1 == METHOD_DCLINDE || method2 == METHOD_DCLINDE) {
    printf("Stage 1 score threshold: %f\n", st1);
    printf("Stage 2 score threshold: %f\n", st2);
    printf("Maximum Parents to condition on in stage 2: %d\n", max_n);
    printf("Pruning: %s\n", pruning_names[pruning]);

    if(is_one_delay)
      printf("Retain only one delay for each link after stage 1.\n");
    if(is_keep_zero_delay)
      printf("To keep links with zero delay after stage 2 and before removing dups.\n");
    if(is_no_dup)
      printf("Retain only one delay for each link, remove duplicate links after stage 2.\n");

    printf("Lower bound (inclusive) of delays to test: %d\n", min_delay);
  }
  printf("Upper bound (inclusive) of delays to test: %d\n", max_delay);

  if(method1 == METHOD_GLOBALMIT || method1 == METHOD_GLOBALMIT_FAST) {
    printf("For initial GRN using GlobalMIT+/*, learn parents starting from gene %d\n", istart);
    if(iend >= 0)
      printf("For initial GRN using GlobalMIT+/*, learn parents for genes up to gene %d\n", iend);
  }
  if(method1 == METHOD_GLOBALMIT || method1 == METHOD_GLOBALMIT_FAST ||
     method2 == METHOD_GLOBALMIT || method2 == METHOD_GLOBALMIT_FAST) {
    printf("Significance level for GlobalMIT+/* (alpha): %g\n", alpha);
    if(is_no_self_link)
      printf("Disallow self links for GlobalMIT+/*.\n");
  }

  if(!is_init_only) {
    printf("Number of states of hidden common cause(s): ");
    if(nh_states < 2)
      printf("Use the max. number of states of children.\n");
    else
      printf("%d\n", nh_states);

    printf("Number of EM iterations: %d\n", n_EM_iter);
    printf("No change iterations for EM restart: %d\n", n_no_change_restart);

    printf("Number of potential parents to find for hidden cause after the first EM: %d\n", h_max_parents);
    printf("Score threshold of association for finding potential parents of hidden cause: %g\n", hp_score_threshold);
  }

  if(s_out_init_grn)
    printf("Output init GRN to file %s\n", s_out_init_grn);
  if(!is_init_only && s_out_grn)
    printf("Output final GRN to file %s\n", s_out_grn);
  /*********/

  data_lens = Malloc(n_data, int, "data_lens");
  data = read_n_data(n_data, argv+s_idx+1, data_lens);
  if(data == NULL) return 1;
  report_stop_watch();

  printf("=== Rename the states if applicable\n");
  rename_data_states(data);
  /* debug 
  {
    printf("data_lens"); pr_ints(n_data, data_lens);
    printf("Concatenated Segments:\n");
    print_carray2d(stdout, data);
  }
  */

  mt_srand(seed);

  /* to prevent abort when there is error in GSL  */
  gsl_set_error_handler(no_abort_handler);
  /*  gsl_set_error_handler_off(); */ 
  /*********/
  hcc_dclinde(data, n_data,data_lens,
	      min_delay, max_delay, nh_states,
	      n_EM_iter, n_no_change_restart,
	      method1, method2,
	      st1,st2, max_n, pruning, 
	      is_one_delay, is_no_dup, is_keep_zero_delay,
	      istart, iend, alpha, is_no_self_link ? 0 : 1,
	      is_init_only, is_no_cpt,
	      eV, eV_tolerance, eV_threshold,
	      h_max_parents, hp_score_threshold,
	      NULL, s_out_grn,s_out_init_grn);

  report_stop_watch();

  /*********/
  free_array2d(data);
  free(data_lens);
  free_temp_buffer();

  if(out_links) free(out_links);

  return 0;
}
