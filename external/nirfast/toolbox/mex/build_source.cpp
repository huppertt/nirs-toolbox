/*
* Modified by: Hamid Ghadyani, March 2009
*           Most of the input arguments are not duplicated within mex file anymore.
*           Fixed seg. fault crashes under Windows 
*/
#include "bem.h"

#define nodes_reg(i,j) nodes_reg[(i)+(nnod_reg*(j))]

void build_source(double *nodes_reg, double *inreg, long nnod, long nnod_reg, mxArray *omega,
double D, double *source, long reg, mxArray *q) {

	mxArray *source_strength, *Gns, *zz;
	mxArray *temp_result;
	double xi, yi, zi; 
	long inod, i;
	double rq;
	
	double *q_real, *q_img;
	double *zz_real, *zz_img, *omega_real, *omega_img, *gns_real, *gns_img;
	double *s_strength_real, *s_strength_img, *temp_real, *temp_img, temp_var;

	source_strength = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
	assert(source_strength != NULL);

	(*(mxGetPr(source_strength))) = 12.4253;
	(*(mxGetPi(source_strength))) = 1.8779;

	zz = mxCreateDoubleMatrix(1,1,mxCOMPLEX);
	Gns = mxCreateDoubleMatrix(1,1,mxCOMPLEX);
	temp_result = mxCreateDoubleMatrix(1,1,mxCOMPLEX);
	
	q_real = mxGetPr(q);
	q_img = mxGetPi(q);

	for(i=0; i<nnod; i++)
		*(q_real+i) = *(q_img+i) = 0.0;	

	zz_real = mxGetPr(zz);
	zz_img = mxGetPi(zz);
	omega_real = mxGetPr(omega);
	omega_img = mxGetPi(omega);

	s_strength_real = mxGetPr(source_strength);
	s_strength_img = mxGetPi(source_strength);
	temp_real = mxGetPr(temp_result);
	temp_img = mxGetPi(temp_result);




	for(i = 0; i < nnod_reg; ++i) {
		inod = (long) inreg[i];

		xi =    nodes_reg(i,0);
		yi =    nodes_reg(i,1);
		zi =    nodes_reg(i,2);
		
		if(reg == 1) {
			rq = sqrt(pow(source[0]-xi,2) + pow(source[1]-yi,2) 
				+ pow(source[2]-zi,2));
			
			*zz_real = (*omega_real) * rq;
			*zz_img = (*omega_img) * rq;
			
			temp_var = 4*PI*rq*D;
			gns_real = mxGetPr(Gns);
			gns_img = mxGetPi(Gns);
			cis((*zz_real * -1.0), (*zz_img * -1.0), *gns_real, *gns_img);
			(*gns_real) /= temp_var;
			(*gns_img) /= temp_var; 

			mult_complex(*gns_real, *gns_img, *s_strength_real, *s_strength_img,
					*temp_real, *temp_img); 
			
			(*(q_real+(inod-1))) += *temp_real;
			(*(q_img+(inod-1))) += *temp_img;
		}

	}
	
		mxDestroyArray(zz); 
		mxDestroyArray(temp_result);
		mxDestroyArray(Gns);
		mxDestroyArray(source_strength);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	long nnod, nnod_reg, reg;
	double *inreg;
	double *source;
	double *nodes_reg;
	double D;
	size_t elt_size;
	mxArray *q; 

	if(nrhs != 7) 
		mexErrMsgTxt("7 inputs required\n");
	
	if(nlhs != 1)
		mexErrMsgTxt("1 output required\n");
	
	nnod = (long) *(mxGetPr(prhs[0]));
	
	elt_size = mxGetElementSize(prhs[1]); assert(elt_size != 0);
	nnod_reg = mxGetM(prhs[1]);
	nodes_reg = mxGetPr(prhs[1]);
	
	if(mxGetN(prhs[2]) != 1) mexErrMsgTxt("Need row vector for inreg\n ");
	inreg = mxGetPr(prhs[2]);
		
	mxArray *omega = mxCreateDoubleMatrix(1,1,mxCOMPLEX); 
	*(mxGetPr(omega)) =  *(mxGetPr(prhs[3]));
	if (mxIsComplex(prhs[3]))
		*(mxGetPi(omega)) =  *(mxGetPi(prhs[3])); 
	else
		*(mxGetPi(omega)) =  0.;
	
	
	D = *(mxGetPr(prhs[4]));

	elt_size = mxGetElementSize(prhs[5]), assert(elt_size != 0);
	source = mxGetPr(prhs[5]);

	reg = (long) *mxGetPr(prhs[6]);
	
	q = mxCreateDoubleMatrix(nnod, 1, mxCOMPLEX);
	
	build_source(nodes_reg, inreg, nnod, nnod_reg, omega, D, source, reg, q);
	
	plhs[0] = q;
	
	mxDestroyArray(omega);
}
