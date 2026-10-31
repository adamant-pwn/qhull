/* Modified by Oleksandr Kulkov with AI assistance, 2026-10-31:
   preserve near-collinear vertices in reentrant 2-d convex hulls.
   Qhull license: COPYING.txt; original source: https://www.qhull.org/ */
/* Private state appended to the library-owned 2-d geometry workspace. */
#ifndef qhDEFrobust2d
#define qhDEFrobust2d
#include "libqhull_r.h"

typedef struct qhRobust2dT {
  boolT enabled;
  realT maxcoord, error;
  coordT initial_vertices[6];
} qhRobust2dT;

/* The first (dim+1)*dim coordinates retain their original scratch role.
   Only the 2-d allocation has a state suffix. Both allocation and release
   use the same size; restarts reuse the buffer and reinitialize this state. */
#define qh_gm_bytes(qh) (((qh)->hull_dim+1)*(qh)->hull_dim*(int)sizeof(coordT) \
    + ((qh)->hull_dim == 2 ? (int)sizeof(qhRobust2dT) : 0))
#define qh_robust2d_state(qh) ((qhRobust2dT *)(void *)((qh)->gm_matrix+6))
#define qh_robust2d(qh) ((qh)->hull_dim == 2 && (qh)->gm_matrix \
    && qh_robust2d_state(qh)->enabled)
#include "exact_sum.h"
/* Resolve exact incidence before deciding whether the original operation
   needs a robust retry. Nonzero near-boundary cases remain conservative. */
static int qh_retry_sign(const coordT *a,const coordT *b,const coordT *c) {
  hull_sum s={{0},{0},0,0};
  hull_product(&s,a[0],b[1],0);hull_product(&s,b[0],c[1],0);hull_product(&s,c[0],a[1],0);
  hull_product(&s,a[1],b[0],1);hull_product(&s,b[1],c[0],1);hull_product(&s,c[1],a[0],1);
  return hull_sum_sign(&s);
}
/* A safe 2-d merge removes only the common point between two collinear
   retained endpoints. Exact collinearity alone permits folded-back edges. */
static boolT qh_local_merge_ok(qhT *qh,facetT *f,facetT *g) {
  vertexT *fv[2], *gv[2];
  coordT *a,*b,*c;
  int i,j,shared=0,fi=0,gi=0;
  if(qh_setsize(qh,f->vertices)!=2 || qh_setsize(qh,g->vertices)!=2)return False;
  for(i=0;i<2;i++) {fv[i]=(vertexT*)SETelem_(f->vertices,i);gv[i]=(vertexT*)SETelem_(g->vertices,i);}
  for(i=0;i<2;i++)for(j=0;j<2;j++)if(fv[i]==gv[j]){shared++;fi=i;gi=j;}
  if(shared!=1)return False;
  a=fv[1-fi]->point;b=fv[fi]->point;c=gv[1-gi]->point;
  if((a[0]==b[0] && a[1]==b[1]) || (c[0]==b[0] && c[1]==b[1]) ||
     (a[0]==c[0] && a[1]==c[1]) || qh_retry_sign(a,b,c))return False;
  for(i=0;i<2;i++)if((b[i]<a[i] && b[i]<c[i]) || (b[i]>a[i] && b[i]>c[i]))return False;
  return True;
}
#endif
