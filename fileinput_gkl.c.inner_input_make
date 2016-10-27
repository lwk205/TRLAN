#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>

#define TWO (2.0)
#define ONE (1.0)
#define ZERO (0.0)

int main(int argc, char *argv[]){

  int s,u,i,each, ii, j, m, n,tmp, w, L ,K, accuracy, lwork;
  double alpha;
  char mode;
  int *IAP,*JA;
  double *A,*work;

  L=atoi(argv[1]);
  accuracy = atoi(argv[3]);
  mode = argv[2][0];
  if(mode=='d') printf("runnning: dense mode\n");
  else if(mode=='s') printf("runnning: sparse mode\n");
  else {printf("error: d[dense] か s[sparse] を指定して下さい。"); return 0;}

  K=L*2;
  scanf("%d",&m);
  scanf("%d",&n);
  scanf("%d",&w);

  printf("each row has %d elements.\n",w/m);

  s=7827187;
  srand(s);
 
  IAP=(int *)malloc(sizeof(int)*(m+1));
  JA=(int *)malloc(sizeof(int)*w);
  A=(double *)malloc(sizeof(double)*w);

  lwork=5*K;
  if (m>lwork) lwork=m;
  if (n>lwork) lwork=n;
  if (K*K>lwork) lwork=K*K;

  work=(double *)malloc(sizeof(double)*lwork);

  if(IAP==NULL || JA==NULL || A==NULL || work==NULL){
    printf("Out of memory.\n");
    return 0;
  }

  each = w / m;
  if(mode=='s'){
    ii=0;
    IAP[0]=1;
    for(tmp=0;tmp<w;tmp++){
    //while(scanf("%d %d %lf",&i,&j,&alpha) !=EOF){
      i = tmp / each;
      j = rand() % n;
      alpha = rand()/(double)RAND_MAX;
      //printf("%d %d %30.20f\n",i,j,alpha);
      if(i!=ii){
        IAP[ii+1]=tmp+1;
        ii++;
      }
      JA[tmp]=j+1;
      A[tmp]=alpha;
      //tmp++;
    }
    IAP[ii+1]=tmp+1;
  }
  else if(mode=='d'){
    tmp=0;
    while(scanf("%lf",&alpha) !=EOF){
      A[tmp]=alpha;
      tmp++;
    }
  }

  resgkl_main_(&mode,&accuracy,&m,&n,&L,&K,IAP,JA,A,work,&lwork);

  return 0;
}

