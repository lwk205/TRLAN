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
  double alpha;
  char mode;
  int *IAP,*JA;
  double *A,*work;

  FILE *fp;
  L=atoi(argv[1]);
  accuracy = atoi(argv[3]);
  mode = argv[2][0];
  fp = fopen(argv[4],"r");
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

  lwork=5*K;
  if (m>lwork) lwork=10*m;
  if (n>lwork) lwork=10*n;
  if (K*K>lwork) lwork=10*K*K;

  work=(double *)malloc(sizeof(double)*lwork);

  if(IAP==NULL || JA==NULL || A==NULL || work==NULL){
    printf("Out of memory.\n");
    return 0;
  }
  
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
  resgkl_main_(&mode,&accuracy,&n,&L,&K,IAP,JA,A,work,&lwork);

  return 0;
}
