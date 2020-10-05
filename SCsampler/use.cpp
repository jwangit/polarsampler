#include "use.h"


void sort(double *data,UINT n,int *index)//sort the data,the array "index" 
{                                                         //store the index of origin data 
	double temp1;
	int temp2;
	for(int j=0;j<n-1;j++)
		for(int i=0;i<n-1-j;i++)
		{
			if(data[i]>data[i+1])//up
			{
			temp1=data[i];
			data[i]=data[i+1];
			data[i+1]=temp1;
			temp2=index[i];
			index[i]=index[i+1];
			index[i+1]=temp2;
			}
		}
}

bool sort(double *data,UINT n,UINT *index)//sort the data,the array "index" 
{                                                         //store the index of origin data 
	double temp1;
	int temp2;
	bool sign=false;
	for(int j=0;j<n-1;j++)
		for(int i=0;i<n-1-j;i++)
		{
			if(data[i]>data[i+1])//up
			{
				sign=true;
			temp1=data[i];
			data[i]=data[i+1];
			data[i+1]=temp1;
			temp2=index[i];
			index[i]=index[i+1];
			index[i+1]=temp2;
			}
		}
	return sign;
}

void sort(double *data,long long n)//sort the data,the array "index" 
{                                                         //store the index of origin data 
	double temp1;
	//int temp2;
	for(long long j=0;j<n-1;j++)
		for(long long i=0;i<n-1-j;i++)
		{
			if(data[i]>data[i+1])//up
			{
			temp1=data[i];
			data[i]=data[i+1];
			data[i+1]=temp1;
			}
		}
}

void sort(UINT *data,UINT n)//sort the data,the array "index" 
{                                                         //store the index of origin data 
	//double temp1;
	int temp2;
	for(UINT j=0;j<n-1;j++)
		for(UINT i=0;i<n-1-j;i++)
		{
			if(data[i]>data[i+1])//up
			{
			temp2=data[i];
			data[i]=data[i+1];
			data[i+1]=temp2;
			}
		}
}

static int graniset = 0;
static double grangset;

double gran(void)
{
   double rsq, v1, v2, fac;

   if(!graniset) {
	  graniset = 1;
	  do {
		 v1 = 2*(rand()/(double)RAND_MAX) - 1;
		 v2 = 2*(rand()/(double)RAND_MAX) - 1;
		 rsq = v1*v1 + v2*v2;
	  } while(rsq > 1 || rsq == 0);  
	  fac = sqrt(-2*log(rsq)/rsq);
	  grangset = v1*fac;
	  return v2*fac;
   }
   else {
	  graniset = 0;
	  return grangset;
   }
}

void randvec(double *y,int N, double sigma)
// fill y with mean= -1 and variance sigma^2 Gaussian
// This corresponds to the all-zero vector
{
   for(int i = 0; i < N; i++) {
	  y[i]=-1+sigma*gran();   //
	
   }
}