#ifndef _USE_H_
#define _USE_H_

#include "stdlib.h"
#include "math.h"
#define PI 3.1415926
typedef unsigned int UINT;

void sort(double *data,UINT n,int *index);//sort the data,the array "index" 

bool sort(double *data,UINT n,UINT *index);//sort the data,the array "index" 

void sort(double *data,long long n);//sort the data,the array "index" 

void sort(UINT *data,UINT n);

double gran(void);

void randvec(double *y,int N, double sigma);
// fill y with mean= -1 and variance sigma^2 Gaussian
// This corresponds to the all-zero vector
#endif