// SCdecoder.cpp 

#include "polar.h"
#include "BECconstruct.h"
#include "time.h"
#include <fstream>
#include <iomanip>
#include <chrono>
#include <mat.h>
#include <matrix.h>
#include <iostream>
#include<sstream> 
#include <fstream>
#include <cstring>
#include <string>
#include <random>
#include <time.h>
#include <numeric>
using namespace std;

std::random_device rd;
std::default_random_engine eng(rd()); // a core engine class
uniform_int_distribution<int> unif(0, 3);

#define Blk 1024//1024 8192 16384  32768: block length N
#define Lev 6		//: number of levels
#define width 64	//: width of a distribution after tailcut, 2^Lev
//#define SCLLR

int main()
{

	double beta = 0.42;		// 2^(-N^beta) is used to define high, lwo, non-polarized sets
	int sigma = 3;			// sigma = s/sqrt(2*pi)
	int n = log2(Blk);
		 
	/*********off-line phase *********/
	/*read channel transition probability table;*/
        /*the transition probabilistic values for level i are stored in transP[i]*/
	const char *InpDistr = "InpDistri_sigmas=3.mat";
	MATFile * pmatFile1 = matOpen(InpDistr, "r");
	if (pmatFile1 == NULL)
	{
		printf("Error opening file %s\n!", InpDistr);
	}
	mxArray *pCell = matGetVariable(pmatFile1, "InpDistri");
	mxArray *pCellArray[Lev];
	double *transP[Lev];
	UINT row = 0; UINT col = 0;
	unsigned int **setflag = new unsigned int*[Lev];
	for (unsigned int i = 0; i < Lev; i++)
	{
		pCellArray[i] = mxGetCell(pCell, i);//mxGetCell_730(pCell, i);
		row = mxGetM(pCellArray[i]);
		col = mxGetN(pCellArray[i]);
		transP[i] = (double*)mxGetPr(pCellArray[i]);

		setflag[i] = new unsigned int[Blk];// remeber to: delete [] setflag[i]  and delete setflag
	}

	matClose(pmatFile1);


	//read Bhattacharyya parameters which ared used to define high, low, non-polarized entropy sets later.
	const string Batafilepre = ("Bata_BDMC_sigmas="+to_string(sigma))+"_Level";
	string level;
	MATFile *pmatFile = NULL;
	mxArray *pMxArray = NULL;
	double *initA;
	for (size_t i = 0; i < Lev; i++)
	{
		level = to_string((i+1));
		string strBatafile = (((Batafilepre + level) + "_n")+to_string(n))+".mat";
		const char * Batafile = strBatafile.c_str();
		cout << strBatafile << endl;
		pmatFile = matOpen(Batafile, "r");//matOpen("Bata_BDMC_sigmas=3_Level4_n10.mat", "r");
		if (pmatFile == NULL) {
			printf("Error opening mat file %s\n", Batafile);
		}
		pMxArray = matGetVariable(pmatFile, "BataLast");
		initA = (double*)mxGetData(pMxArray);
		row = mxGetM(pMxArray);
		col = mxGetN(pMxArray);
		printf(matClose(pmatFile) == 0 ? "mat file closed\n" : "Error: mat file closed\n");

		double delta = pow(2, -pow(Blk, beta));
		constructPolar(initA, delta, Blk, setflag[i]);		 // define high low and unpolarized sets using delta.


/*			mxDestroyArray(pMxArray);
			mxFree(pMxArray);
			mxDestroyArray(pCell);
			mxFree(pCellArray[0]);
		    mxDestroyArray(pCellArray[0]);*/

	}

	double *W0 = new double[Blk];
	double *W1 = new double[Blk];
	double *Wm0 = new double[Blk];
	double *Wm1 = new double[Blk];
	UINT * X = new UINT[Blk];

	UINT N = Blk;
	UINT M = N / 2;
	UINT m = (UINT)log2(N);
	double numiter = 10000;     //number of iterations;to invoke polar sampler repeatedly
	UINT * preX = new UINT [Blk] ;
	UINT **itransP = new UINT*[Lev];

	//Since the values in LR[i] and transP[i] are stored in bitreversal order of X0X1...X(i-1),
        // we setup the one-to-one mapping between itransP[i][j] and X0X1...X(i-1). 
	for (UINT i = 0; i < Lev; i++)
	{
		itransP[i] = new UINT [int(pow(2, i))];
			for (UINT j = 0; j < int(pow(2, i)); j++)
				itransP[i][j] = j;
		bitReversal(itransP[i], int(pow(2, i)), i);	
	}
	srand(time(NULL));
	polar polarAll[Lev] = {
	polar(N, M, m, setflag[0]),
	polar(N, M, m, setflag[1]),
	polar(N, M, m, setflag[2]),
	polar(N, M, m, setflag[3]),
	polar(N, M, m, setflag[4]),
	polar(N, M, m, setflag[5]),
	};
	UINT histArray[width] = { 0 };
	
	//Since high entropy set is uniformly at random, it can be produced offline.
	polarAll[0].randhighEnt();
	polarAll[1].randhighEnt();
	polarAll[2].randhighEnt();
	polarAll[3].randhighEnt();
	polarAll[4].randhighEnt();
	polarAll[5].randhighEnt();
	double speed = 0;
	
	/********online phase********/
	auto start = chrono::steady_clock::now();
	for (UINT J = 0; J < numiter; J++)
	{
		for (size_t i = 0; i < Blk; i++)
		{
			preX[i] = 0;// unif(eng);
		}
		for (UINT k = 0; k < Lev; k++)
		{
			for (UINT i = 0; i < Blk; i++)
			{

				W0[i] = *(transP[k] + 2 * itransP[k][preX[i]]);
				W1[i] = *(transP[k] + 2 * itransP[k][preX[i]] + 1);

			}

			//polarAll[k].randhighEnt();

			polarAll[k].SCdecoder(W0, W1, Wm0, Wm1, eng);
			X = polarAll[k].encode();
			for (UINT i = 0; i < Blk; i++)
			{
				preX[i] += (X[i] * pow(2, k));
			}

		}
		/******statistics**********/
		/*for (UINT i = 0; i < Blk; i++)
		{
			histArray[preX[i]] ++;
		}*/
	}
	/*******speed output *******/
	auto end = chrono::steady_clock::now();
	speed = Blk*numiter / chrono::duration_cast<chrono::microseconds>(end - start).count() * 1000000;
	cout << "elapsed time in us:" << chrono::duration_cast<chrono::microseconds>(end - start).count() << " us" << endl;
	cout.precision(15);
	cout << "speed:" << speed << " samples/s" << endl;
	/******statistics output**********/
	/*cout << endl;
	cout << "Histogram value at index..." << endl;
	for (UINT i = 0; i < pow(2,Lev); i++)
	{
		cout << i << ":" << histArray[i] / (Blk * numiter)<< endl;
	}*/

	

	return 0;
}
