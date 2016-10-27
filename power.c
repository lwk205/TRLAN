#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
double gettimeofday_sec(){
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec+(double)tv.tv_usec*1e-6;
}
#define TWO (2.0)  
#define ONE (1.0)
#define ZERO (0.0)
int M, N;
void sparse_matvec(int *IAP,int *JA,double *A,double *p,
                   double *Ap){
  int i,j;
  for(i=0;i<M;i++){
    Ap[i]=ZERO;
    for(j=IAP[i];j<IAP[i+1];j++){
      Ap[i]=Ap[i]+(A[j])*p[JA[j]];
    }
  }
  return;
}
void sparse_matvec2(int *IAP,int *JA,double *A,double *q,
                    double *Aq){
  int i,j;
  for(i=0;i<N;i++){
    Aq[i]=ZERO;
  }
  for(i=0;i<M;i++){
    for(j=IAP[i];j<IAP[i+1];j++){
      Aq[JA[j]]=Aq[JA[j]]+(A[j])*q[i];
    }
  }
  return;
}
int main(int argc, char *argv[]){
  int i, j, k, l;
  int W, MPN;
  double alpha, norm_x, norm_r, sigma;
  double t1, t2;
  srand(1);
  W=atoi(argv[1]);
  scanf("%d",&M);
  scanf("%d",&N);
  scanf("%d",&l);
  MPN=M+N;
  double x[W][N], y[W][M], Ax[M], ATy[N], r[MPN];
  int IAP[M+1], JA[l];
  double A[l];
  
  k = 0;
  l = 0;
  IAP[0] = 0;
  while(scanf("%d %d %lf",&i,&j,&alpha) !=EOF){
    if(i!=k){
      IAP[++k] = l;
    }
    JA[l] = j;
    A[l] = alpha;
    l++;
  }
  IAP[k+1] = l;

  t1=gettimeofday_sec();
  for(i=0;i<W;i++){
    for(k=0;k<N;k++){ x[i][k]=TWO*rand()/(RAND_MAX)-ONE; }
    for(j=0;j<100*MPN;j++){
      /* modified Gram-Schmidt */
      for(k=0;k<i;k++){
        norm_r=ZERO;
        for(l=0;l<N;l++){ norm_r += x[k][l]*x[i][l]; }
        for(l=0;l<N;l++){ x[i][l] -= norm_r*x[k][l]; }
      }
      /* norm = 1 */
      norm_x=ZERO;
      for(k=0;k<N;k++){ norm_x += x[i][k]*x[i][k]; }
      norm_x=sqrt(norm_x);
      for(k=0;k<N;k++){ x[i][k] /= norm_x; }
      sparse_matvec(IAP,JA,A,&(x[i][0]),Ax);
      for(k=0;k<M;k++){ y[i][k]=Ax[k]; }
      /* modified Gram-Schmidt */
      for(k=0;k<i;k++){
        norm_r=ZERO;
        for(l=0;l<M;l++){ norm_r += y[k][l]*y[i][l]; }
        for(l=0;l<M;l++){ y[i][l] -= norm_r*y[k][l]; }
      }
      /* norm = 1 */
      norm_x=ZERO;
      for(k=0;k<M;k++){ norm_x += y[i][k]*y[i][k]; }
      norm_x=sqrt(norm_x);
      for(k=0;k<M;k++){ y[i][k] /= norm_x; }
      sparse_matvec2(IAP,JA,A,&(y[i][0]),ATy);
      /* stopping criterion */
      sigma=ZERO;
      for(k=0;k<N;k++){ sigma += x[i][k]*ATy[k];  }
      printf("j=%d sigma=%f\n",j,sigma);
      for(k=0;k<M;k++){ r[k]=Ax[k]-sigma*y[i][k]; }
      for(k=0;k<N;k++){ r[M+k]=ATy[k]-sigma*x[i][k]; }
      /* residual error */
      norm_r=ZERO;
      for(k=0;k<MPN;k++){ norm_r += r[k]*r[k]; }
      norm_r=sqrt(norm_r/TWO);
      if(norm_r<1.0E-8){ break; }
      for(k=0;k<N;k++){ x[i][k] = ATy[k]; }
    }
    if (j==100*MPN){ printf("not converge\n"); return 0; }
    printf("sigma=%f\n",sigma);
    fflush(stdout);
  }
  t2=gettimeofday_sec();
  printf("TIME=%10.5g\n",t2-t1);
  return 0;
}
