#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

/* #include "mt_rand.h" */
#include "parse_option.h"

/**************************************************************
  To accept a weighted (undirected) graph, and to detect possibly
  overlapping communities, by means of finding semi-cliques and patch
  up the rest.

***********************************************************************/
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

/****************************************************************/
#define LINK_FROM   1
#define LINK_TO     2
#define LINK_DELAY  4
#define LINK_EFFECT 8
/****************************************************************/
typedef struct s_edge {
  int from, to; /* 0-based indices */
  int delay;
  double effect;
} edge;

/* For temporary use in read_grn() */
static int n_edges=0;
static int n_vertices=0; /* 1 + the maximum index seen in the non-ignored edges */
static int n_capacity_edges=0;
static edge* all_edges=NULL;

static void add_edge(int from, int to, int delay, double effect) {
  if(n_edges >= n_capacity_edges) { /* expand buffer */
    n_capacity_edges = n_capacity_edges < 10 ? 10 : 2*n_capacity_edges;
    all_edges = (edge*)realloc(all_edges, n_capacity_edges*sizeof(edge));
    if(all_edges == NULL) {
      fprintf(stderr, "Not enough memory for %d edges in add_edge().\n", n_capacity_edges);
      exit(-1);
    }
  }

  all_edges[n_edges].from = from;
  all_edges[n_edges].to = to;
  all_edges[n_edges].delay = delay;
  all_edges[n_edges].effect = effect;
  n_edges++;

  if(from+1 > n_vertices) n_vertices = from+1;
  if(to+1 > n_vertices) n_vertices = to+1;
}

static void clear_all_edges() {
  if(all_edges) free(all_edges);
  all_edges = NULL;
  n_capacity_edges = 0;
  n_edges = 0;
  n_vertices = 0;
}

static int read_word(FILE* f, char* buf, int L) {
  /* read next item and until next whitespace, read at most L-1
     characters, and put '\0' at the end.

     return EOF if end of file reached before reading any thing
     interesting.
     return 0 if end of line reached before reading other things.
     Otherwise returns the number of characters read.
  */
  int i=0, c=0;
  while((c=fgetc(f))==' ' || c=='\t' || c=='\r')
    ;
  if(c==EOF) return EOF;
  if(c=='\n') return 0;

  buf[i] = c;
  for(i=1; i<L-1; i++) {
    c = fgetc(f);
    if(c==' ' || c=='\t' || c=='\r' || c==EOF)
      break;
    if(c=='\n') {ungetc(c,f); break;}
    buf[i] = c;
  }
  buf[i] = '\0';

  /*
  printf("==== debug: read_word(): \"%s\"\n", buf);
  */

  return i;
}
static int read_one_edge(FILE* f, int* from, int* to, int* delay, double* effect) {
  /* either
       From: from To: to Delay: delay Coef: effect
     or
       From: from To: to Delay: delay Effect: effect

     The order of from, to, delay and coef (effect)
     can be different, other things on the same line are ignored.

     return EOF if reached end of file before reading anything
     interesting.  otherwise returns the union of flags of LINK_FROM,
     LINK_TO, LINK_DELAY, LINK_EFFECT, to indicates which are
     successfully read.
 */
  int r=0, flag=0;
#define MAX_BUF 64
  static char buf[MAX_BUF]; /* should be large enough for our use */
  r = read_word(f, buf, MAX_BUF);
  if(r == EOF) return EOF;
  if(r == 0) return 0; /* read nothing in this line */
  /* keyword form */
  while(1) {
    for(r=0; buf[r]!='\0'; r++) buf[r] = tolower(buf[r]);
    if(0==strcmp(buf, "from:")) {
      if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
      *from = atoi(buf);
      flag |= LINK_FROM;
    } else if(0==strcmp(buf, "to:")) {
      if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
      *to = atoi(buf);
      flag |= LINK_TO;
    } else if(0==strcmp(buf, "delay:")) {
      if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
      *delay = atoi(buf);
      flag |= LINK_DELAY;
    } else if(0==strcmp(buf, "coef:") || 0==strcmp(buf, "effect:")) {
      if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
      *effect = atof(buf);
      flag |= LINK_EFFECT;
    } else { /* ignore other keywords and its associated value */
      if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
    }
    /* possibly next keyword */
    if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
  }

  return flag;
}
static void read_edges_from_file(char* name, int required) {
  FILE* f=NULL;
  int from=0, to=0, delay=0, n=0, i=0;
  double effect=0;

  f = fopen(name, "r");
  if(f == NULL) {
    fprintf(stderr, "Fail to read \"%s\" for edges.\n", name);
    exit(-1);
  }

  i=0;
  while(!feof(f)) {
    delay=0; effect=0;
    n = read_one_edge(f, &from, &to, &delay, &effect);
    if(n==EOF || (n == 0 && feof(f))) {
      break;
    } else if(n==0) {
      continue;
    } else if((n & required) != required) {
      /* some required field missing */
      fprintf(stderr, "Error in reading edge file, after reading %d edge(s).\nIn the edge file, each line represents an edge, with the 0-based 'from' index, followed by 0-based 'to' index, followed by delay, followed by the effect, all separated by space.\n", i);
      if((n & LINK_FROM)==0)   fprintf(stderr,"Link missing 'from'\n");
      if((n & LINK_TO)==0)     fprintf(stderr,"Link missing 'to'\n");
      if((n & LINK_DELAY)==0)  fprintf(stderr,"Link missing 'delay'\n");
      if((n & LINK_EFFECT)==0) fprintf(stderr,"Link missing 'effect' or 'coef'\n");
      exit(-1);
    }
    add_edge(from, to, delay, effect);
    i++;
  }

  fclose(f);
}

/* A simple wrapper */
edge* read_grn(char* name, int* out_n_edges, int* out_n_vertices, int* out_max_dely, int required) {
  /* To read the edges and gives some useful information. The read
     edges should be freed after use. */
  int i=0, d=0;
  edge* r=NULL;
  clear_all_edges();
  read_edges_from_file(name, required);
  if(out_n_edges) *out_n_edges = n_edges;
  if(out_n_vertices) *out_n_vertices = n_vertices;

  if(out_max_dely) {
    for(i=0; i<n_edges; i++) {
      if(all_edges[i].delay > d) d = all_edges[i].delay;
    }
    *out_max_dely = d;
  }
  /* transfer ownership of all_edges to here. */
  r = all_edges;
  all_edges = NULL;
  clear_all_edges();

  return r;
}

static void print_grn(FILE* f, int n_edges, edge* es) {
  /* mainly for debug */
  int i=0;
  fprintf(f, "=========================================\n");
  fprintf(f, "%d edges:\n", n_edges);
  for(i=0; i<n_edges; i++)
    fprintf(f, "%d:\tfrom: %d\tto: %d\tdelay: %d\teffect: %f\n", i,
	    es[i].from, es[i].to, es[i].delay, es[i].effect);
  fprintf(f, "=========================================\n");
}

/***************************************************************/
int ascending_from_idx(const void* a, const void* b) {
  /* *a and *b contain edge* */
  edge *pa=NULL, *pb=NULL;
  int r=0;
  pa = (*(edge**)a);
  pb = (*(edge**)b);
  r = pa->from - pb->from;
  if(r != 0) return r;
  return pa->delay - pb->delay;
}

int ascending_to_idx(const void* a, const void* b) {
  /* *a and *b contain edge* */
  edge *pa=NULL, *pb=NULL;
  int r=0;
  pa = (*(edge**)a);
  pb = (*(edge**)b);
  r = pa->to - pb->to;
  if(r != 0) return r;
  return pa->delay - pb->delay;
}
/***************************************************************/
typedef struct s_vertex {
  int n_from;
  edge** from; /* not own this, pointer to the edges with the same from */
  int n_to;
  edge** to;  /* not own this, pointer to the edges with the same to */
} vertex;

void index_grn(int n_edges, edge* es, int n_vertices, vertex* vs, edge* buf[]) {
  /* will write to vs. buf should have length 2*n_edges, and is used
     as buffer in vs. */
  int i=0, k=0;
  edge** p=NULL;

  for(i=0; i<n_vertices; i++) {
    vs[i].n_from = 0;
    vs[i].n_to = 0;
    vs[i].from = NULL;
    vs[i].to = NULL;
  }

  /* count the number of edges for each vertex first */
  for(i=0; i<n_edges; i++) {
    vs[es[i].from].n_from++;
    vs[es[i].to].n_to++;
  }

  /* distribute the buffers */
  for(i=0, p=buf; i<n_vertices; i++) {
    vs[i].from = p;
    p += vs[i].n_from;
    vs[i].n_from = 0; /* to use as index to count again */

    vs[i].to = p;
    p += vs[i].n_to;
    vs[i].n_to = 0; /* to use as index to count again */
  }

  /* debug 
  printf("**** n_edges: %d\tp-buf: %d\n", n_edges, p-buf);
  fflush(stdout);
  */

  /* point appropriately */
  for(i=0; i<n_edges; i++) {
    k = es[i].from;
    vs[k].from[vs[k].n_from++] = es+i;

    k = es[i].to;
    vs[k].to[vs[k].n_to++] = es+i;
  }

  /* sort the links */
  for(i=0; i<n_vertices; i++) {
    qsort(vs[i].from, vs[i].n_from, sizeof(edge*), ascending_to_idx);
    qsort(vs[i].to, vs[i].n_to, sizeof(edge*), ascending_from_idx);
  }
}

void print_indexed_grn(FILE* f, int n_vertices, vertex vs[]) {
  int i=0, k=0;
  edge** p = NULL;

  fprintf(f,"**** index GRN\n");
  for(i=0; i<n_vertices; i++) {
    fprintf(f,"vertex %d (%d + %d edges):\n", i, vs[i].n_from, vs[i].n_to);
    p = vs[i].from;
    for(k=0; k<vs[i].n_from; k++) {
      fprintf(f,"\t=> %d delay %d effect %f\n", p[k]->to, p[k]->delay, p[k]->effect);
    }
    p = vs[i].to;
    for(k=0; k<vs[i].n_to; k++) {
      fprintf(f,"\t<= %d delay %d effect %f\n", p[k]->from, p[k]->delay, p[k]->effect);
    }
  }
  fprintf(f,"****\n");
}

/***************************************************************/
static int sign_of(double x) {
  if(x < 0) return -1;
  if(x > 0) return 1;
  return 0;
}

static double cal_fmeasure(double r, double p) {
  if(r < 1.0e-8 || p < 1.0e-8) return 0;
  return 2*r*p/(r+p);
}

static void n_from_correct(int n1, edge* es1[], int n2, edge* es2[],
		     int* nc_links, int* nc_delay, int* nc_effect) {
  /* how many edges of es1 that are same as in es2, for links, delay
     and effect.

     The edges in es1 and es2 are assumed to be from the same
     vertex, sorted by ascending to index (then by ascending delay).

     Output to nc_links, nc_delay, nc_effect.
   */
  int nL=0, nD=0, nE=0;
  int i=0, j=0, k=0, delay_ok=0, effect_ok=0, s=0;
  edge *e1=NULL, *e2=NULL;

  /* for debug
  printf("***** %d compare with %d edges\n", n1, n2);
  for(i=0; i<n1; i++)
    printf("\tes1: => %d, d %d, %f\n", es1[i]->to,es1[i]->delay,es1[i]->effect);
  for(j=0; j<n2; j++)
    printf("\tes2: => %d, d %d, %f\n", es2[j]->to,es2[j]->delay,es2[j]->effect);
  printf("***\n");
  */

  i=0; j=0;
  while(i<n1 && j<n2) {
    e1 = es1[i];
    e2 = es2[j];

    /* for debug
    printf("\tes1[%d]%d,d %d,%f\tes2[%d]%d,d %d,%f\n",
	   i,e1->to,e1->delay,e1->effect,
	   j,e2->to,e2->delay,e2->effect);
    */

    if(e1->to < e2->to) {
      i++;
    } else if(e1->to == e2->to) {
      nL++;
      /* there may be multiple edges for the same to, use brute
	 force for delay and effect */
      delay_ok = 0;
      effect_ok = 0;
      s = sign_of(e1->effect);
      for(k=j; k<n2 && es2[k]->to == e1->to; k++) {
	if(es2[k]->delay == e1->delay) delay_ok = 1;
	if(sign_of(es2[k]->effect) == s) effect_ok = 1;
      }
      nD += delay_ok;
      nE += effect_ok;
      i++; /* j would be incremented when the e1->from changes */
    } else { /* e1->to > e2->to */
      j++;
    }
  }

  if(nc_links) *nc_links = nL;
  if(nc_delay) *nc_delay = nD;
  if(nc_effect) *nc_effect = nE;
}

void align_node_to_node(vertex* a, vertex* b, int n_non_hidden,
			int max_delay, int buf[],
			int* score, int* d_shift, int* flip) {
  /* align a to b, return the best delay match possible as score.
     d_shift: amount to increase for parent links, to decrease for
              children links, could be negative.
     flip: whether to flip the sign of the effects of all links

     The output score is the sum of number of matched links and
     delays.  buf should have length 1+2*max_delay: the first
     max_delay will be for the negative shifts, then 0 shift, then
     positive shifts.
  */
  int i=0, j=0, k=0;
  int nL=0, d=0, s=0, s_matched=0, s_not_matched=0,
    n_e_matched=0, n_e_not_matched=0;
  edge *ea=NULL, *eb=NULL, *e=NULL;

  for(i=0; i<(1+2*max_delay); i++) buf[i] = 0; /* for voting which delay shift has most matches */

  /* only look at those neighbors < n_non_hidden, because the hidden
     nodes may have been renamed.

     Mainly walk through a.
  */
  /* for 'from', children */
  i=0; j=0;
  while(i<(a->n_from) && j<(b->n_from)) {
    ea = a->from[i];
    eb = b->from[j];
    if(ea->to >= n_non_hidden || eb->to >= n_non_hidden) break;

    if((ea->to) < (eb->to)) {
      i++;
    } else if((ea->to) == (eb->to)) {
      nL++;
      s = sign_of(ea->effect);
      s_matched = 0; s_not_matched = 0;
      for(k=j, e=b->from[k]; k<(b->n_from) && e->to == ea->to; k++, e=b->from[k]) {
	if(sign_of(e->effect) == s) s_matched=1; else s_not_matched=1;
	d = ea->delay - e->delay; /* shift in children is to be subtracted */
	buf[max_delay + d]++;
      }
      n_e_matched += s_matched; n_e_not_matched += s_not_matched;
      i++; /* j will be incremented when ea->to changes */
    } else { /* (ea->to) > (eb->to) */
      j++;
    }
  }
  /* similarly for 'to', parents */
  i=0; j=0;
  while(i<(a->n_to) && j<(b->n_to)) {
    ea = a->to[i];
    eb = b->to[j];
    if(ea->from >= n_non_hidden || eb->from >= n_non_hidden) break;

    if((ea->from) < (eb->from)) {
      i++;
    } else if((ea->from) == (eb->from)) {
      nL++;
      s = sign_of(ea->effect);
      s_matched = 0; s_not_matched = 0;
      for(k=j, e=b->to[k]; k<(b->n_to) && e->from == ea->from; k++, e=b->to[k]) {
	if(sign_of(e->effect) == s) s_matched=1; else s_not_matched=1;
	d = e->delay - ea->delay; /* shift in parent is to be added */
	buf[max_delay + d]++;
      }
      n_e_matched += s_matched; n_e_not_matched += s_not_matched;
      i++; /* j will be incremented when ea->from changes */
    } else { /* (ea->from) > (eb->from) */
      j++;
    }
  }

  /* now determine the best delay shift and effect sign flip */
  for(i=0, j=0, d=0; i<(1+2*max_delay); i++) {
    if(buf[i] > d) {
      d = buf[i]; /* best number of matched delay so far */
      j = i;
    }
  }
  if(d_shift) *d_shift = j - max_delay;

  if(flip) {
    *flip = (n_e_not_matched > n_e_matched && flip) ? 1 : 0;
  }

  if(score) *score = nL + d;
}

void align_hidden_node(vertex* h, int n, vertex vs[], int n_non_hidden,
		       int max_delay, int buf[]) {
  /* possibly adjust the to, from, delay and sign of effect in h.
     buf should have length 1+2*max_delay, for temporary use.
   */
  int best_idx=0, best_score=0, best_d_shift=0, best_flip=0;
  int score=0, d_shift=0, flip=0;
  int i=0;
  edge* e=NULL;

  best_idx = n_non_hidden;
  for(i=n_non_hidden; i<n; i++) {
    score=0; d_shift=0; flip=0;
    align_node_to_node(h,vs+i, n_non_hidden,max_delay,buf,
		       &score,&d_shift,&flip);

    if(score > best_score) {
      best_idx = i;
      best_score = score;
      best_d_shift = d_shift;
      best_flip = flip;
    }
  }
  /* modify h to align */
  /* d_shift: amount to increase for parent links, to decrease for
              children links, could be negative.
     flip: whether to flip the sign of the effects of all links
  */

  /* for debug 
  printf("\taligned to true node %d, score: %d, d_shift: %d, flip: %d\n",
	 best_idx, best_score, best_d_shift, best_flip);
  */

  for(i=0; i<(h->n_to); i++) { /* parent */
    e = h->to[i];
    e->to = best_idx;
    e->delay += best_d_shift;
    if(best_flip) e->effect = -(e->effect);
  }
  for(i=0; i<(h->n_from); i++) { /* children */
    e = h->from[i];
    e->from = best_idx;
    e->delay -= best_d_shift;
    if(best_flip) e->effect = -(e->effect);
  }
}

void cmp_grn(char* pred_f, char* true_f, int n_non_hidden, int verbose, int required) {
  /* prints the links.recall, links.precision, links.fmeasure,
     delay.recall, delay.precision, delay.fmeasure,
     effect.recall, effect.precision, effect.fmeasure.

     But required is a union of flags indicating the required fields
     among from, to, delay and effect: LINK_FROM, LINK_TO, LINK_DELAY
     and LINK_EFFECT.

     If n_non_hidden <= 0, means no hidden nodes. Otherwise all
     indices >= n_non_hidden are regarded as hidden nodes.
  */
  int i=0, n=0, max_delay=0;
  int nL_pred=0, nv_pred=0, md_pred=0;
  int nL_true=0, nv_true=0, md_true=0;
  double nL_r=0, nD_r=0, nE_r=0, nL_p=0, nD_p=0, nE_p=0;
  int nL=0, nD=0, nE=0;
  edge *pred_grn=NULL, *true_grn=NULL;
  vertex *pred_vs=NULL, *true_vs=NULL;
  edge** buf=NULL;
  int* shift_buf=NULL;

  /* read predicted grn */
  pred_grn = read_grn(pred_f, &nL_pred, &nv_pred, &md_pred, required);
  if(verbose) {
    printf("Predicted GRN, %d edges, %d genes, max-delay: %d.\n",
	   nL_pred, nv_pred, md_pred);
    print_grn(stdout, nL_pred, pred_grn);
    printf("\n");
  }
  /* read true grn */
  true_grn = read_grn(true_f, &nL_true, &nv_true, &md_true, required);
  if(verbose) {
    printf("True GRN, %d edges, %d genes, max-delay: %d.\n",
	   nL_true, nv_true, md_true);
    print_grn(stdout, nL_true, true_grn);
    printf("\n");
  }

  max_delay = md_pred > md_true ? md_pred : md_true;
  /* to prepare at least nv_true slots for pred_vs, because after aligning hidden nodes,
     the number of vertices of pred may increase. */
  n = nv_pred > nv_true ? nv_pred : nv_true;

  /* index them for more efficient processing */
  buf = Malloc(2*(nL_pred+nL_true), edge*, "cmp_grn()");
  memset(buf, 0, sizeof(edge*)*2*(nL_pred+nL_true));
  shift_buf = Malloc(1+2*max_delay, int, "cmp_grn()");
  memset(shift_buf, 0, sizeof(int)*(1+2*max_delay));
  pred_vs = Malloc(n, vertex, "cmp_grn()");
  memset(pred_vs, 0, sizeof(vertex)*n);
  true_vs = Malloc(nv_true, vertex, "cmp_grn()");
  memset(true_vs, 0, sizeof(vertex)*nv_true);

  index_grn(nL_pred, pred_grn, nv_pred, pred_vs, buf);
  index_grn(nL_true, true_grn, nv_true, true_vs, buf+2*nL_pred);

  if(verbose) {
    printf("**** Predicted GRN\n");
    print_indexed_grn(stdout, nv_pred, pred_vs);
    printf("**** True GRN\n");
    print_indexed_grn(stdout, nv_true, true_vs);
    printf("\n");
  }
  /* debug 
  fflush(stdout);
  */

  /* align the predicted hidden nodes to the true hidden nodes */
  if(n_non_hidden > 0 && nv_true > n_non_hidden && nv_pred > n_non_hidden) {
    for(i=n_non_hidden; i<nv_pred; i++) {
      /* for debug 
      printf("Aligning predicted node %d\n", i); fflush(stdout);
      */
      align_hidden_node(pred_vs+i, nv_true,true_vs,
			n_non_hidden, max_delay, shift_buf);
    }

    /* re-index the predicted GRN, to merge hidden nodes that are
       aligned to the same true hidden node. */
    index_grn(nL_pred, pred_grn, n, pred_vs, buf); /* number of vertices may have increased */
    if(verbose) {
      printf("**** Predicted GRN after aligning hidden nodes\n");
      print_indexed_grn(stdout, n, pred_vs);
      printf("\n");
    }
  }

  /* get recall, precision, need only look at from */
  if(nv_true < n) n = nv_true;
  for(i=0; i<n; i++) {
    /* for debug
    printf("\nvertex %d recall\n", i);
    */
    /* recall */
    nL=0; nD=0; nE=0;
    n_from_correct(true_vs[i].n_from, true_vs[i].from,
		   pred_vs[i].n_from, pred_vs[i].from,
		   &nL, &nD, &nE);
    nL_r += nL; nD_r += nD; nE_r += nE;

    /* for debug
    printf("\nvertex %d precision\n", i);
    */
    /* precision */
    nL=0; nD=0; nE=0;
    n_from_correct(pred_vs[i].n_from, pred_vs[i].from,
		   true_vs[i].n_from, true_vs[i].from,
		   &nL, &nD, &nE);
    nL_p += nL; nD_p += nD; nE_p += nE;
  }
  /* calculate F-measure */
  if(nL_true > 0) {
    nL_r /= nL_true; nD_r /= nL_true; nE_r /= nL_true;
  } else {
    nL_r = 0; nD_r = 0; nE_r = 0;
  }
  if(nL_pred > 0) {
    nL_p /= nL_pred; nD_p /= nL_pred; nE_p /= nL_pred;
  } else {
    nL_p = 0; nD_p = 0; nE_p = 0;
  }
  if(verbose) {
    printf("==== Comparison\n==== link.r,link.p,link.f");
    if(required & LINK_DELAY)
      printf(",delay.r,delay.p,delay.f");
    if(required & LINK_EFFECT)
      printf(",effect.r,effect.p,effect.f");
    printf("\n");
  }
  /* links */
  printf("%f,%f,%f", nL_r,nL_p,cal_fmeasure(nL_r,nL_p));
  if(required & LINK_DELAY)
    printf(",%f,%f,%f",	nD_r,nD_p,cal_fmeasure(nD_r,nD_p));
  if(required & LINK_EFFECT)
    printf(",%f,%f,%f",	nE_r,nE_p,cal_fmeasure(nE_r,nE_p));
  printf("\n");

  /* cleanup */
  free(pred_grn);
  free(true_grn);
  free(buf);
  free(shift_buf);
  free(pred_vs);
  free(true_vs);
}

/***************************************************************/
#define DEFAULT_NON_HIDDEN      0

#define xstr(s) str(s)
#define str(s) #s

/* The option structure of this program */
option all_options[] = {
  /* str,type(STRING),is_required(0),description,index(0),end_index(0), val_int(0),val_real(0),val_str(NULL) */
  {"-?",      FLAG,   0,  "Showing the usage.\n", 0,0,  0,  0.0,    NULL},

  {"-p",   STRING,    1,  "File name of the predicted GRN. Each line in the file represents an edge, which consists of a 0-based 'from' index, a 0-based 'to' index, a delay, and the effect, separated by space..\n", 0,0,  -1,  0.0,    NULL},
  {"-t",   STRING,    1,  "File name of the true GRN. Each line in the file represents an edge, which consists of a 0-based 'from' index, a 0-based 'to' index, a delay, and the effect, separated by space..\n", 0,0,  -1,  0.0,    NULL},
  {"-n",   INT,       0,  "Number of non-hidden genes, unspecified if <= 0. If an index is >= this, it is regarded as a hidden node. Default " xstr(DEFAULT_NON_HIDDEN) ".\n", 0,0,  DEFAULT_NON_HIDDEN,  0.0,    NULL},
  {"-no.effect",    FLAG,   0,  "No effect.\n", 0,0,  0,  0.0,    NULL},
  {"-no.delay",     FLAG,   0,  "No delay.\n", 0,0,  0,  0.0,    NULL},
  {"-v",      FLAG,   0,  "Verbose mode.\n", 0,0,  0,  0.0,    NULL},

  {NULL,      STRING, 0,  NULL, 0,0,  0,  0.0,    NULL} /* end marker */
};
/***************************************************************/

int main(int argc, char *argv[])
{
  time_t starting_time, end_time;
  char *pred_grn_file_name=NULL, *true_grn_file_name=NULL;
  int n_non_hidden = DEFAULT_NON_HIDDEN;
  int verbose=0, no_effect=0, no_delay=0;
  int required=0;

  required = LINK_FROM | LINK_TO | LINK_DELAY | LINK_EFFECT;

  /* do some initialization */
  if(parse_options(argc,argv,all_options)) {
    /* error parsing options */
    usage(stderr, argv[0], all_options);
    exit(0);
  }
  if(get_flag_option_value("-?",all_options,0))
    usage(stdout, argv[0], all_options);
  verbose = get_flag_option_value("-v",all_options,0);
  no_effect = get_flag_option_value("-no.effect",all_options, no_effect);
  no_delay = get_flag_option_value("-no.delay",all_options, no_delay);

  if(no_effect) required &= (~LINK_EFFECT);
  if(no_delay)  required &= (~LINK_DELAY);

  /*********************/
  if(verbose) {
    time(&starting_time);
    printf("starting time: %ld\t%s", (long)starting_time, ctime(&starting_time));
    if(no_effect) printf("Ignore effects in the links.\n");
    if(no_delay) printf("Ignore delays in the links.\n");
  }

  /*********************/
  n_non_hidden = get_int_option_value("-n",all_options, n_non_hidden);
  if(verbose) printf("Number of non-hidden nodes: %d\n", n_non_hidden);

  pred_grn_file_name = get_str_option_value("-p",all_options,NULL);
  if(verbose) printf("Predicted GRN file name: %s\n", pred_grn_file_name);
  true_grn_file_name = get_str_option_value("-t",all_options,NULL);
  if(verbose) printf("True GRN file name: %s\n", true_grn_file_name);

  cmp_grn(pred_grn_file_name, true_grn_file_name, n_non_hidden, verbose, required);

  /*********************/

  if(verbose) {
    time(&end_time);
    printf("ending time: %d\t%s", (long)end_time, ctime(&end_time));
    printf("Total seconds used: %d\n", end_time - starting_time);
    printf("====End\n\n");
  }

  return 0;
}
