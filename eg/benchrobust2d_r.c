/* Oleksandr Kulkov with AI assistance, 2026-10-31.
   Measure qh_new_qhull plus cleanup, excluding input parsing.
   Qhull license: COPYING.txt. Original source: https://www.qhull.org/ */
#include "libqhull_r/qhull_ra.h"
#include <time.h>
int main(int argc,char **argv) {
  int dim,n,i,j,code,c,l,vertices=0,repeat=argc>1?atoi(argv[1]):5;
  coordT *points;
  clock_t start=0,elapsed=0;
  FILE *quiet=tmpfile();
  if(!quiet || repeat<1 || scanf("%d%d",&dim,&n)!=2 || dim!=2 || n<3)return 1;
  points=(coordT*)malloc((size_t)n*2*sizeof(coordT));if(!points)return 1;
  for(i=0;i<2*n;i++)if(scanf("%lf",points+i)!=1)return 1;
  for(j=-1;j<repeat;j++) {
    qhT state,*qh=&state;
    if(j==0)start=clock();
    qh_zero(qh,quiet);
    code=qh_new_qhull(qh,dim,n,points,False,"qhull Qt",NULL,quiet);
    if(code)return code;
    vertices=qh->num_vertices;
    qh_freeqhull(qh,!qh_ALL);qh_memfreeshort(qh,&c,&l);
    if(c||l)return 2;
  }
  elapsed=clock()-start;
  printf("%.9f %d\n",1000.*elapsed/CLOCKS_PER_SEC/repeat,vertices);
  free(points);fclose(quiet);return 0;
}
