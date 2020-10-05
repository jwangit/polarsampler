#include "polar.h"
#include <iomanip>
#include<fstream>
#include <random>
#include <time.h>

ofstream write("u.txt");
polar::polar(UINT N,UINT M,UINT m, UINT *setflag_)
{
	this->N=N;
	this->M=M;
	this->m=m;
	K=N-M;
	R=double(K)/N;
	alloc();
/*	for(UINT i=0;i<M;++i)
		frozen[i]=frozen_[i];
	for(UINT i=0;i<N;++i)
		frozenflag[i]=0;
	for(UINT i=0;i<M;++i)
		frozenflag[frozen[i]]=1;*/
	for (UINT i = 0; i < N; i++)
		setflag[i] = setflag_[i];
}
polar::~polar()
{
	free();
}
void polar::alloc()
{
	u=new UINT [N];
	ub=new UINT [N];
	x=new int [N];
	y=new double [N];
	//bitB=new UINT [N];
	frozen=new UINT [M];
	frozenflag=new UINT [N];
	setflag = new UINT[N];
	arrayPointer_LR = new double *[m + 1];
	arrayPointer_P=new double **[m+1];
	arrayPointer_C=new UINT **[m+1]; 
	for(UINT i=0;i<=m;++i)
    { 
		UINT size=(1<<(m-i));
	    arrayPointer_P[i]=new double *[size];
	    arrayPointer_C[i]=new UINT *[size];
		arrayPointer_LR[i] = new double[size];
	    for(UINT j=0;j<size; ++j)
		{ 
			arrayPointer_P[i][j]=new double[2];// dynamically allocate transition probability array and bits array
		    arrayPointer_C[i][j]=new UINT[2];
		}
    }
}
void polar::free()
{
	delete [] u;
	delete [] ub;
	delete [] x;
	delete [] y;
//	delete [] bitB;
	delete [] frozen;
	delete [] frozenflag;
	delete[] setflag;
	for(UINT i=0;i<=m;++i)
	{
		UINT size=(1<<(m-i));
	    for(UINT j=0;j<size;++j)
		{ 
			delete [] arrayPointer_P[i][j];
		    delete [] arrayPointer_C[i][j];
		}
		delete [] arrayPointer_P[i];
	    delete [] arrayPointer_C[i];
		delete[] arrayPointer_LR[i];
	}
	delete [] arrayPointer_P;
	delete [] arrayPointer_C;
	delete[] arrayPointer_LR;
}
void polar::randomVec()
{
	cout<<"u sequence:"<<endl;
	write<<"u sequence:"<<endl;
	for(unsigned int i=0;i<N;++i)
	{
		if(frozenflag[i]==1)
			u[i]=0;
		else
			
			u[i]=rand()%2;
		cout<<u[i]<<"  ";
		write<<"u["<<i<<"]="<<u[i]<<"  ";
		if((i+1)%8==0)
		{
			cout<<endl;	write<<endl;
		}
	}
	cout<<endl;
	write<<endl;write<<endl;write<<endl;write<<endl;
}
void polar::randhighEnt()
{

	for (unsigned int i = 0; i < N; ++i)
	{
		if (setflag[i] == 1)
			u[i] = rand() % 2;
		else
			u[i] = 0;
	}

}
void polar::Reversal() //rearrange ub in bitreversal order e.g.: 1234  -> 1324
{
	UINT *index=new UINT [N];
	for(UINT i=0;i<N;++i)
		index[i]=i;
	bitReversal(index);
	for(UINT i=0;i<N;++i)
		ub[index[i]]=u[i];
/*	cout<<"bitreversal后ub:"<<endl;;
	for(UINT i=0;i<N;++i)
	{
		cout<<ub[i]<<"  ";
		if((i+1)%8==0)
			cout<<endl;
	}
	cout<<endl;*/
	delete [] index;
}
void polar::bitReversal(UINT *index)
{
	UINT *index1=new UINT [N];
	UINT m1=m-1;
	for(UINT i=0;i<m1;++i)
	{
		UINT count1=(1<<i);
		UINT count2=N>>(i+1);
		for(UINT k=0;k<count1;++k)
		{
			UINT len=k*(N>>i);        
			for(UINT j=0;j<count2;++j)
			{
				index1[j+len]=index[2*j+len];
				index1[j+count2+len]=index[2*j+1+len];
			}
		}
		for(UINT j=0;j<=N-1;++j)
			index[j]=index1[j];
	}
	delete [] index1;
}
UINT* polar::encode()
{
	Reversal();
	bool *sign=new bool [N];
	for(UINT i=1;i<=m;++i)
	{
		for(UINT k=0;k<N;++k)
			sign[k]=false;
		for(UINT j=0;j<N;++j)
		{
			if(sign[j])
				continue;
			else
			{
				UINT pos=j+(N>>i);
				ub[j]=ub[j]^ub[pos];
				sign[pos]=true;
			}
		}
	}

	/*cout<<"encoded ub:"<<endl;;
	for(UINT i=0;i<N;++i)
	{
		cout<<ub[i]<<"  ";
		if((i+1)%8==0)
			cout<<endl;
	}
	cout<<endl;*/
	delete [] sign;
	return ub;
}
void polar::BPSKmodulation()
{
	//cout<<"BPSK调制：x值"<<endl;
	for(unsigned int i=0;i<N;++i)
	{
		if(ub[i]==0)
			x[i]=-1;
		else 
			x[i]=1;
	//	cout<<x[i]<<"  ";
		//if((i+1)%8==0)
			//cout<<endl;
	}
}
void polar::channelOut(double sigma) //
{
	//cout<<endl;
	//cout<<"y:"<<endl;
	for(UINT i=0;i<N;++i)
	{
		if(x[i]==1)
			y[i]=1+sigma*gran();
		else
			y[i]=-1+sigma*gran();
	//	cout<<y[i]<<"  ";
		//if((i+1)%4==0)
		//	cout<<endl;
	}
}
void  polar::transp(double *W0,double *W1,double sigma2)  //W(y|x)
{
	double factor = 1.0/(2*PI*sqrt(sigma2));
	double index = -1.0/(2*sigma2);
	for(UINT i=0; i<N;++i)
	{  
		W0[i]=factor*exp(index*pow(y[i]+1,2.0)); 
		W1[i]=factor*exp(index*pow(y[i]-1,2.0));
    }
	
}

void polar::recursivelyCalcLR(UINT lamda, UINT phi)
{
	if (lamda == 0) return;
	UINT pus = phi / 2;

	if (phi % 2 == 0)
		recursivelyCalcLR(lamda - 1, pus);
	double deta = 0;
	UINT **C;
	double *LR1 = arrayPointer_LR[lamda];
	double *LR0 = arrayPointer_LR[lamda - 1];


	UINT size = (1 << (m - lamda));
	for (UINT beta = 0; beta < size; ++beta)
	{
		if (phi % 2 == 0)
		{
			LR1[beta] = (LR0[2 * beta] * LR0[2 * beta + 1] + 1) / (LR0[2 * beta] + LR0[2 * beta + 1]);
			//P1[beta] = log((exp(P0[2 * beta] + P0[2 * beta + 1]) + 1) / (exp(P0[2 * beta]) + exp(P0[2 * beta + 1])));
			LR1[beta] = (LR1[beta] > LRmax) ? LRmax : LR1[beta];
			LR1[beta] = (LR1[beta] < (1 / LRmax)) ? (1 / LRmax) : LR1[beta];
			//deta = (deta >= LR1[beta] ? deta : LR1[beta]);
		}
		else
		{
			C = arrayPointer_C[lamda];
			UINT u1 = C[beta][0];
			LR1[beta] = pow(LR0[2 * beta], (1 - 2 * u1)) * LR0[2 * beta + 1];
			//P1[beta] = (1 - 2 * up) * P1[2 * beta] + P1[2 * beta + 1];
			LR1[beta] = (LR1[beta] > LRmax) ? LRmax : LR1[beta];
			LR1[beta] = (LR1[beta] < (1 / LRmax)) ? (1 / LRmax) : LR1[beta];
		}
	}

}

void polar::recursivelyCalcP(UINT lamda,UINT phi)
{
	 if(lamda==0) return ;
	 UINT pus=phi/2;

	 if(phi%2==0)   
		 recursivelyCalcP(lamda-1,pus);
	 double deta=0;
	 UINT **C;
	 double **P1=arrayPointer_P[lamda];
	 double **P0=arrayPointer_P[lamda-1];
     UINT size=(1<<(m-lamda));
	 for(UINT beta=0;beta<size;++beta)
	 {
		 if(phi%2==0)
			 for(UINT u=0;u<2;++u)
			 {
				 P1[beta][u]=0.5*(P0[2*beta][u]*P0[2*beta+1][0]+P0[2*beta][1^u]*P0[2*beta+1][1]);
				 deta=(deta>=P1[beta][u]?deta:P1[beta][u]);
/*				 if (P1[beta][u] > 1.0 || P1[beta][u] < 0.0)
				 {
					 printf("Error: invalid prob!");
				 }*/
			 }
		 else
		 {
			 C=arrayPointer_C[lamda];
			 UINT u1=C[beta][0];
			 for(UINT u2=0;u2<2;++u2)
			 {
				 P1[beta][u2]=0.5*P0[2*beta][u1^u2]*P0[2*beta+1][u2];
				 deta=(deta>=P1[beta][u2]?deta:P1[beta][u2]);
/*				 if (P1[beta][u2] > 1.0 || P1[beta][u2] < 0.0)
				 {
					 printf("Error: invalid prob!");
				 }*/
			 }
		 }
	 }
	 for(UINT beta=0;beta<size;++beta)
		 for(UINT u=0;u<2;++u)
			 P1[beta][u]=P1[beta][u]/deta;
}
void polar::recursivelyUpdataC(UINT lamda,UINT phi)
{
	if(lamda==0) return;

	UINT pus=phi/2;
	UINT size=1<<(m-lamda);
	UINT **C0=arrayPointer_C[lamda-1];
	UINT **C1=arrayPointer_C[lamda];
	for(UINT beta=0;beta<size;++beta)
	{ 
		C0[2*beta][pus%2]=C1[beta][0]^C1[beta][1];
	    C0[2*beta+1][pus%2]=C1[beta][1];
	}
	if(pus%2==1)	
	   recursivelyUpdataC(lamda-1,pus);
}
void polar::SCdecoder(double *W0,double *W1,double *Wm0,double *Wm1, std::default_random_engine eng)
{
	double **P0=arrayPointer_P[0];
	double **Pm=arrayPointer_P[m];
	UINT **Cm=arrayPointer_C[m];
    for(UINT beta=0;beta<N;++beta)
	{
		P0[beta][0]=W0[beta];
	    P0[beta][1]=W1[beta];
	}
/*	for(UINT phi=0;phi<N;++phi)
	{
		recursivelyCalcP(m,phi);
		Wm0[phi]=Pm[0][0];
		Wm1[phi]=Pm[0][1];

		if(frozenflag[phi]==1)
		{
			Cm[0][phi % 2] = u[phi];//  //rand() % 2
			//bitB[phi]= u[phi];
		}
		else if(Wm0[phi]>Wm1[phi])
		{
			Cm[0][phi%2]=0;
			//bitB[phi]=0;
		}
		else
		{
			Cm[0][phi%2]=1;
			//bitB[phi]=1;
		}
		if(phi%2==1)
			recursivelyUpdataC(m,phi);
	}*/
	for (UINT i = 0; i < N; ++i)
	{
		recursivelyCalcP(m, i);
		Wm0[i] = Pm[0][0];
		Wm1[i] = Pm[0][1];
		if (setflag[i] == 1)
		{
			Cm[0][i % 2] = u[i];// rand() % 2;
			//bitB[i] = u[i];// rand() % 2;
		}
		else if (setflag[i] == 2)
		{
			Cm[0][i % 2] = Wm0[i] > Wm1[i] ? 0 : 1;
			//bitB[i] = Cm[0][i % 2];
			u[i] = Cm[0][i % 2];
		}
		else
		{
//			std::random_device rand_ber;
//			std::default_random_engine engine(rand_ber()); // a core engine class
			double temp1 = Wm1[i] / (Wm0[i] + Wm1[i]);
			std::bernoulli_distribution bern(temp1);
			Cm[0][i % 2] = bern(eng) ? 1 : 0 ;
			u[i] = Cm[0][i % 2];
			//bitB[i] = Cm[0][i % 2];
		}
		if (i % 2 == 1)
			recursivelyUpdataC(m,i);
	}
	//return u;
	//	return bitB;
}
void polar::SCdecoder_LR(double * LR, double * LRm, std::default_random_engine eng)
{
	double *Plr0 = arrayPointer_LR[0];
	double *Plrm = arrayPointer_LR[m];
	UINT **Cm = arrayPointer_C[m];
	for (UINT beta = 0; beta < N; ++beta)
	{
		Plr0[beta] = LR[beta];
	}
	for (UINT i = 0; i < N; ++i)
	{
		recursivelyCalcLR(m, i);
		LRm[i] = Plrm[0];

		if (setflag[i] == 1)
		{
			Cm[0][i % 2] = u[i];// rand() % 2;
			//bitB[i] = u[i];// rand() % 2;
		}
		else if (setflag[i] == 2)
		{
			Cm[0][i % 2] = (LRm[i] > 1) ? 0 : 1;
			//bitB[i] = Cm[0][i % 2];
			u[i] = Cm[0][i % 2];
		}
		else
		{
			double temp = 1 / (LRm[i] + 1);
			if (temp < 0.0 || temp > 1.0)
			{
				printf("\n Error: invalid Bernoulli parameter!\n");
			}
			std::bernoulli_distribution bern(temp);
			Cm[0][i % 2] = bern(eng) ? 1 : 0;
			u[i] = Cm[0][i % 2];
			
		}
		if (i % 2 == 1)
			recursivelyUpdataC(m, i);

	}
}
//void polar::callDecoder(double *W0,double *W1,double *Wm0,double *Wm1)
//{
//	SCdecoder(W0,W1,Wm0,Wm1);
//}
/*UINT polar::bitCount()
{
	UINT errCount=0;
	for(UINT i=0;i<N;++i)
		if(frozenflag[i]==1||bitB[i]==u[i])
			continue;
		else
			++errCount;
	//cout<<"errCount="<<errCount<<endl;
	return errCount;
}*/

/*void polar::statistics(UINT count0, UINT count1, UINT numiter)
{
	for (UINT i = 0; i < N; i++)
	{
		if (ub[i] == 0)
		{
			count0 += 1;
		}
		else if (ub[i] == 1)
		{
			count1 += 1;
		}
		else
		{
			cout << "SCdecoder error: ub[i] invalid. \n" << endl;
			break;
		}
	}

}*/


