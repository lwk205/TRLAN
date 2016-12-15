#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>

#define TWO (2.0)
#define ONE (1.0)
#define ZERO (0.0)

int main(int argc, char *argv[]){

  int i, ii, j, m, n, w, L ,K, accuracy, lwork;
  int maxthreads,k0,k1,t;
  double alpha,n0,l0,l1;
  char mode;
  int *IAP,*JA,*start_row;
  double *A,*work;
  int which;
  FILE *fp;
  which = atoi(argv[1]);
  // 1 == 絶対値最大からK
  // 2 == 絶対値最小からK
  // 3 == 値最大からK
  // 4 == 値最小からK
  L=atoi(argv[2]);
  mode = argv[3][0];
  accuracy = atoi(argv[4]);
  fp = fopen(argv[5],"r");
  if(mode=='d') printf("runnning: dense mode\n");
  else if(mode=='s') printf("runnning: sparse mode\n");
  else {printf("error: d[dense] か s[sparse] を指定して下さい。"); return 0;}

  K=L*2;
  fscanf(fp,"%d",&m);
  fscanf(fp,"%d",&n);
  fscanf(fp,"%d",&w);

  printf("each row has %d elements.\n",w/m);

  IAP=(int *)malloc(sizeof(int)*(m+1));
  JA=(int *)malloc(sizeof(int)*w);
  A=(double *)malloc(sizeof(double)*w);

  lwork=5*K*20;
  if (m>lwork) lwork=m;
  if (n>lwork) lwork=n;

  work=(double *)malloc(sizeof(double)*lwork);
  
  maxthreads=omp_get_max_threads();
  start_row=(int *)malloc((maxthreads+1)*sizeof(int));

  if(IAP==NULL || JA==NULL || A==NULL || work==NULL || start_row==NULL){
    printf("Out of memory.\n");
    return 0;
  }
  
  n0 = (double)w/maxthreads;

  if(mode=='s'){
    ii=0;
    w=0;
    IAP[0]=1;
    while(fscanf(fp,"%d %d %lf",&i,&j,&alpha) !=EOF){
      if(i!=ii){
        IAP[ii+1]=w+1;
        ii++;
      }
      JA[w]=j+1;
      A[w]=alpha;
      w++;
    }
    IAP[ii+1]=w+1;
  }
  else if(mode=='d'){
    w=0;
    while(fscanf(fp,"%lf",&alpha) !=EOF){
      A[w]=alpha;
      w++;
    }
  }
  fclose(fp);


  for(t=0,i=0;i<maxthreads && t<m;i++){
    start_row[i]=t;

    k1=0;
    while(k1<n0 && t<m){
      k0=k1;
      k1=k0+IAP[t+1]-IAP[t];
      t=t+1;
    }

    if(start_row[i]+1==t){
      continue;
    }

    l0=n0-k0;
    l1=fabs(k1-n0);
    if(l0<l1){
      t=t-1;
    }
  }
  start_row[i]=m;

  if(i!=maxthreads){
    for(i=i+1;i<maxthreads+1;i++){
      start_row[i]=0;
    }
  }
  for(i = 0;i<maxthreads+1;i++){
  printf("start_row %d = %d\n",i,start_row[i]);
  }

  resl_main_(start_row,&which,&mode,&accuracy,&m,&L,&K,IAP,JA,A,work,&lwork);

  return 0;
}
