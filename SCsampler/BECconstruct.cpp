#include "BECconstruct.h"


void bitReversal(unsigned int *ivec,int len,int pow)
{
	int pow2=pow;//log2(len);
	
	int *tmpvec=new int [pow2];
	int tmpint=0;
	unsigned mask=1;
	int i,j;
	
	for (i=0;i<len;i++)
	{
		tmpint=ivec[i];
		mask=1;
		mask<<=pow2;
		for (j=0;j<pow2;j++)
		{
			tmpint<<=1;
			tmpvec[j]=mask&tmpint;
		}
		
		mask=1;
		for (j=0;j<pow2;j++)
		{
			if (tmpvec[j])
			{
				ivec[i] |= mask;
			} 
			else
			{
				ivec[i] &= (~mask);
			}
			mask<<=1;
		}
	}
	delete [] tmpvec;
}




void constructPolar(double *initA, double delta, unsigned int iN, unsigned int *setflag)//, unsigned int *frozen, unsigned int *lowEntro)//setup high and low entropy index flag
{
	unsigned int N = iN;
	int i = 0;
	unsigned int *index = new unsigned int[N];
	double *Zn = new double[N];

	for (unsigned int j = 0; j < N; j++) 
	{
		Zn[j] = initA[j];
		index[j] = j;
	}
	sort(Zn, N, index);
	for (unsigned int i = 0; i < N; i++)
	{	
		if (Zn[i] < delta)
		{
			setflag[index[i]] = 2;

		}
		else if (Zn[i] > 1-delta)
		{
			setflag[index[i]] = 1;
		}
		else
		{
			setflag[index[i]] = 0;
		}
	}

	delete[] Zn;
	delete[] index;
}
