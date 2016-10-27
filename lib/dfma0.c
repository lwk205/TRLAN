#include<stdio.h>
#include<math.h>
#define FP_FAST_FMA

double dfma0_(double *x,double *y,double *z){
  return fma(*x,*y,*z);
}
