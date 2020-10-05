#include<iostream>
#include "use.h"
#include <random>
using namespace std;

class polar
{
public:
	polar(UINT N,UINT M,UINT m, UINT *setflag_);
	~polar();
	void alloc();
	void free();
public:
	void randomVec();
	void randhighEnt();
	void Reversal();
	void bitReversal(UINT *index);
	//void encode();
	UINT * encode();
	void BPSKmodulation();
	void channelOut(double sigma);
	void transp(double *W0,double *W1,double sigma2);  //W(y|x)
private:
	void recursivelyCalcLR(UINT lamda, UINT phi);
	void recursivelyCalcP(UINT lamda,UINT phi);
	void recursivelyUpdataC(UINT lamda,UINT phi);
//	void SCdecoder(double *W0,double *W1,double *Wm0,double *Wm1);
public:
	void SCdecoder(double *W0, double *W1, double *Wm0, double *Wm1, std::default_random_engine  eng);
	void SCdecoder_LR(double *LR, double *LRm, std::default_random_engine eng);
//	void callDecoder(double *W0,double *W1,double *Wm0,double *Wm1);
//	UINT bitCount();
//	void statistics(UINT count0, UINT count1, UINT numiter);

private:
	UINT N,M,K;
	UINT m;
	double R;
	double LRmax=1E30;
private:
	UINT *u,*ub;
	UINT bitB;
	int *x;
	double *y;
private:
	UINT *frozen;      //frozen position
	UINT *frozenflag;
	UINT *setflag;
private:
	 double ***arrayPointer_P;
     unsigned int ***arrayPointer_C;
	 double **arrayPointer_LR;
};

// std::random_device rd; 
// std::mt19937 eng(rd());  // a core engine class
//extern std::tr1::uniform_int<int> unif0(0, 1);
//extern std::tr1::uniform_int<int> unif1(0, 3);

//std::tr1::bernoulli_distribution bern(0.25);