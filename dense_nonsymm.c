#include<stdio.h>
#include<stdlib.h>
int main(int argc,char *argv[]){
  int i, j, l, m, n, w;
  m=atoi(argv[1]);
  n=atoi(argv[2]);
  srand(1);
  double a;
  printf("%d\n",m);
  printf("%d\n",n);
  printf("%d\n",n*m);
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      a=rand()/(double)RAND_MAX;
      printf("%30.20f\n",a);
    }
  }
  return 0;
}
