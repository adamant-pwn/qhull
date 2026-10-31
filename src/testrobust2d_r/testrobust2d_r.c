/* Modified by Oleksandr Kulkov with AI assistance, 2026-10-31:
   regression tests for near-collinear two-dimensional hull vertices.
   Distributed under Qhull's license in COPYING.txt; original source:
   https://www.qhull.org/ */
#include "libqhull_r/qhull_ra.h"
#include <math.h>

static int check(const double *input, int n, int expected, const char *label) {
  qhT state, *qh=&state;
  coordT points[8];
  vertexT *vertex;
  int code, curlong, totlong, count=0, mask=0, i, failed=0;
  FILE *quiet=tmpfile();
  if (!quiet) return 1;
  for (i=0;i<2*n;i++) points[i]=input[i];
  qh_zero(qh,quiet);
  code=qh_new_qhull(qh,2,n,points,False,"qhull Qt",NULL,quiet);
  if (!code) {
    FORALLvertices {
      int id=qh_pointid(qh,vertex->point);
      if (id<0 || id>=n) failed=1;
      else mask |= 1<<id;
      count++;
    }
  }
  if (expected) failed |= code || count!=expected || mask!=((1<<n)-1);
  else failed |= !code;
  qh_freeqhull(qh,!qh_ALL);
  qh_memfreeshort(qh,&curlong,&totlong);
  failed |= curlong!=0 || totlong!=0;
  fclose(quiet);
  if (failed) fprintf(stderr,"%s: exit=%d vertices=%d expected=%d\n",label,code,count,expected);
  return failed;
}
int main(void) {
  const double base[8]={0,0,100000000,99999999,-100000000,100000000,99999999,99999998};
  const double line[8]={0,0,1,1,2,2,3,3};
  const double box[8]={0,0,3,0,3,4,0,4};
  double p[8];
  int a,b,c,d,i,e,failed=0;
  for(a=0;a<4;a++)for(b=0;b<4;b++)for(c=0;c<4;c++)for(d=0;d<4;d++) {
    int order[4];
    if(a==b || a==c || a==d || b==c || b==d || c==d)continue;
    order[0]=a;order[1]=b;order[2]=c;order[3]=d;
    for(i=0;i<4;i++){p[2*i]=base[2*order[i]];p[2*i+1]=base[2*order[i]+1];}
    failed |= check(p,4,4,"near-collinear permutation");
  }
  for(e=-900;e<=900;e+=900){
    for(i=0;i<8;i++)p[i]=ldexp(box[i],e);
    failed |= check(p,4,4,"scaled rectangle");
  }
  failed |= check(line,4,0,"collinear input");
  return failed;
}
