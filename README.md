# polarsampler
A discrete Gaussian sampler over the integers 

1. A polar sampler consists of an online phase and an offline phase. The offline phase is finished in MATLAB by which the data for online phase are produced and saved in .mat format. This repo gives a demo of the online phase in C++ language. It focuses on how the online phase works as well as speed tests. The data type of most probabilistic values is made to be double. The C++ library <random> is used to generate basic random numbers. So there is no constant-time guarantee. 

2. File list in folder SCsampler
a. “Bata_BDMC_sigmas=3_Levelx_ny.mat”: Battacharyya parameters computed offline in MATLAB and stored in mat file; x means the i th level of binary  partition; y     indicates block length N=2^y.
b. “InpDistri_sigmas=3.mat”: contains the transition probabilities of each level
c. “BECconstruct.cpp”, “BECconstruct.h ”: define the high and low entropy sets; define the bit-reversal function.
d. “decoding.cpp”: main function
e. “polar.cpp”, “polar.h”: define polar sampler class
f. “use.h”, “use.cpp”: define other functions

3. Commands:
command 1 : $export LD_LIBRARY_PATH=/usr/local/MATLAB/R2016b/bin/glnxa64
command 2 : $g++ -Ofast decoding.cpp BECconstruct.cpp polar.cpp use.cpp -o SCsampler -I/usr/local/MATLAB/R2016b/extern/include -L/usr/local/MATLAB/R2016b/bin/glnxa64 -lmat -leng -lmx -lmex
command 3 : $./SCsampler
To build and run the codes, you are suggested to install MATLAB on your PC (any version after MATLAB2016 is supposed to work correctly). Please find the folder "glnxa64" and the folder "include" in the installation directory and use command 1 and 2 to build. Please make sure the path in command 1 and 2 are correct according to your installation directory. 
The reason we use MATLAB library is that it provides interfaces between .mat data and C++. The calculations in offline phase are done in MATLAB and the data including Bhattacharyya parameters and transition probabilities are saved in the MATLAB data format. It allows C++ to read .mat data correctly.
       
