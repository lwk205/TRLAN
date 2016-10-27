#include <stdio.h>
#include <stdlib.h>

#define ZERO (double)0.0
#define TWO (double)2.0

int main(int argc,char *argv[])
{
  long long i, j, k, n, m, z;
  double s;
  n=atoi(argv[1]);
  m=atoi(argv[2])/2;
  srand(1);
  double a[n][n];
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      a[i][j]=ZERO;
    }
  }
  for(i=0;i<n*m;i++){
    j=rand() % n;
    k=rand() % n;
    s=(double)rand()/(RAND_MAX);
    while(a[j][k]!=ZERO){
        j=rand() % n;
        k=rand() % n;
    }
    a[j][k]=s;
    a[k][j]=s;
  }
  z=0;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(a[i][j]!=ZERO)
        z=z+1;
    }
  }
  printf("%lld\n",n);
  printf("%lld\n",n);
  printf("%lld\n",z);// z はだいたいargv[2]と一致
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(a[i][j]!=ZERO){
        printf("%lld %lld %30.20f\n",i,j,a[i][j]);
      }
    }
  }
  return 0;
}
