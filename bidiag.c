#include<stdio.h>
#include<stdlib.h>
int main(int argc,char *argv[]){
  int i, j, k, l, m, n, w;
  FILE *fp;
  fp = fopen(argv[2],"w");
  m=atoi(argv[1]);
  srand(1);
  double a;
  fprintf(fp,"%d\n",m);
  fprintf(fp,"%d\n",m);
  fprintf(fp,"%d\n",2*m-1);
  for(i=0;i<m-1;i++){
      a=rand()/(double)RAND_MAX;
      fprintf(fp,"%d %d %30.20f\n",i,i,a);
      a=rand()/(double)RAND_MAX;
      fprintf(fp,"%d %d %30.20f\n",i,i+1,a);
  }
      a=rand()/(double)RAND_MAX;
      fprintf(fp,"%d %d %30.20f\n",m-1,m-1,a);
  
  fclose(fp);
  return 0;
}
