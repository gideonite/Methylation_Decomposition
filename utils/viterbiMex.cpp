#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <mex.h>
#include <iostream>
#define pi 3.1415926535

/* inner product of bits a with w */
double inner(int a, double* w, int d){
	double cuml = 0;

	for(int i = 0; i < d; i++)
		cuml = cuml + ((a>>i)&1)*w[i];

	return cuml;
}

/* normal pdf at v */
double normpdf(double v, double mu, double sigma){
	return 1/(sigma*sqrt(2*pi)) * exp(-pow(v-mu,2.0)/(2.0*pow(sigma,2.0)));
}

/* integer power function */
int myPow(int a, int b){
	return (int)pow((double) a, (double) b);
}

/* max of two numbers */
double max(double a, double b){
	return (a>b)?a:b;
}

/* max index of double array */
int maxIndDouble(double l[],int size){
	int max_ind = 0;

	for(int i = 0; i < size; i++){
		if(l[i] > l[max_ind])
			max_ind = i;
	}
	return max_ind;
}

/* max index of double array */
int maxIndInt(int l[],int size){
	int max_ind = 0;
	for(int i = 0; i < size; i++){
		if(l[i] > l[max_ind])
			max_ind = i;
	}
	return max_ind;
}


/* */
int maxMinInd(double* l, int size, double dist){
	for(int i = 0; i < size; i++){
		if(l[i] >= dist){
			return i;		
		}
	}
	return 1;
}


/* log of binomial pdf */
double logbinompdf(double k, double n, double p){
	double sum1 = 0, sum2 = 0;
	for(int i = n-k+1; i <= n; i++)
		sum1 = sum1 + log(i);
	for(int i = 1; i <=k; i++)
		sum2 = sum2 + log(i);

	return sum1 - sum2 + k*log(p) + (n-k)*log(1-p);
}

/*
  TODO do whatever it takes to implement the beta distribution at v,
  just as in the other functions. Will probably involve getting the
  gamma distribution. See
  http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
  for the chance to take some code??
*/

void viterbi_unique(double* x0, double** Ts, double* bins, double* pos, double* w, 
	double noise, int depth, int d, int t, int n, double* s_star){
	int D = myPow(2,d);
	/*double V[n][D];*/
	double** V = new double*[n];
	for(int i = 0; i < n; i++)
		V[i] = new double[D];

	/* initialize to -inf */
	for(int z = 0; z < n; z++){
		for(int q = 0; q < D; q++){
			V[z][q] = -DBL_MAX;
		}
	}

	/* fill in first row */
	for(int i = 0; i < D; i++){
		double stateTw = inner(i,w,d);
		if(noise == 0){
			V[0][i] = d*log(0.5) + logbinompdf(depth*x0[0],depth,stateTw);
		}
		else{
			double sigma = max(sqrt(depth*stateTw*(1-stateTw)*noise), 1e-3);
			V[0][i] = d*log(0.5) + log(normpdf(depth*x0[0],depth*stateTw,sigma));
		}
	}
	
	for(int i = 1; i < n; i++){
		double dist = pos[i] - pos[i-1];
		double* T;
		if(dist <= 0){
			T = (double *)malloc(4*sizeof(double));
			T[0] = 0.5;
			T[1] = 0.5;
			T[2] = 0.5;
			T[3] = 0.5;
			if (dist == 0)
			{
				printf("DISTANCE IS 0 AT %d\n", i);
			}
		}
		else{
			int maxind = maxMinInd(bins,t,dist);
			T = Ts[maxind-1]; /*?????*/
		}
		for(int j = 0; j < D; j++){
			for(int prev = 0; prev < D ; prev++){
				double log_em_prob = 0;
				double stateTw = inner(j,w,d);
				if(noise == 0){
					log_em_prob = logbinompdf(depth*x0[i], depth, stateTw);
				}
				else{
					double sigma = max(sqrt(depth*stateTw*(1-stateTw)*noise), 1e-3);
					log_em_prob = log(normpdf(depth*x0[i], depth*stateTw, sigma));
				}

				double log_tr_prob = 0.0;

				for(int k = 0; k < d; k++){
					int prev_bin = ((prev>>k)&1);
					int j_bin = ((j>>k)&1);
					log_tr_prob = log_tr_prob +log(T[prev_bin + j_bin*2]);
				}

				V[i][j] = max(V[i][j],log_em_prob + log_tr_prob + V[i-1][prev]);
			}
		}
		if(dist < 0){
			free(T);
		}
	}

	for(int j = 0; j < n; j++){
		s_star[j] = maxIndDouble(V[j],D) + 1;
	}
}

/* params
 * 0 x0 	length n double*
 * 1 Ts 	length t cell matrix 
 * 2 bins	length t integer*
 * 3 pos	length n integer*
 * 4 w		length d double*
 * 5 noise	scalar double
 * 6 depth	scalar integer 
 * 7 d		scalar 
 * 8 t		scalar 
 * 9 n		scalar 
 */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){
	
	int d = mxGetScalar(prhs[7]);
	int t = mxGetScalar(prhs[8]);
	int n = mxGetScalar(prhs[9]);

	double* x0 = (double*) mxGetPr(prhs[0]);
	const mxArray* cell = prhs[1];
	double* output = mxGetPr(cell);
	const mwSize* dims;
	mwIndex jcell;
	/*dims = mxGetDimensions(prhs[1]);*/
	double** Ts = (double**)mxMalloc(sizeof(double*)*t);
	for(jcell = 0; jcell < t-1; jcell++){
		Ts[jcell] = (double*)mxGetPr(mxGetCell(cell,jcell));
	} 

	double* bins = (double*) mxGetPr(prhs[2]);
	double* pos = (double*) mxGetPr(prhs[3]);
  	double* w = (double*) mxGetPr(prhs[4]);
	double noise = mxGetScalar(prhs[5]);
	double depth = mxGetScalar(prhs[6]);

	plhs[0] = mxCreateDoubleMatrix(1,n,mxREAL);
	double* S_star = mxGetPr(plhs[0]);
	viterbi_unique(x0,Ts, bins, pos, w, noise, depth, d, t, n, S_star);
	
}
