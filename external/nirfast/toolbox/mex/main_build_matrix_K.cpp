/**
* main_build_matrix.c
* FUNCTIONS: 
*	main_build_matrix 	-Entry point of console application
*
* Serves as entry point for application
* 
* Written By: 
*	Senate Taka from a translation of a fortran-90 program
* 	that does the same thing. Original fortran code by Subha
*
* For:
*	NIR Dartmouth College 2008
*
* Created: 05 Sept 2008
* Modified: Hamid Ghadyani
*				March 2010
*           	  Most of the input arguments are not duplicated within mex file anymore.
*           	  Fixed seg. fault crashes under Windows
*				Apr 2010: 
*				  Cleaned up usage of macros to access matrices.
*				  Removed unnecessary copying of elements of integral results.
*				May 2010:
*				  Added OpenMP to significantly reduce the runtime
*/



/* How to mex
*	if your compiler supports OpenMP, use this:
*   GCC:
*   mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" main_build_matrix_K.cpp
*
*	Intel/{Mac,Linux}
*   mex -v CXXFLAGS="\$CXXFLAGS -openmp" LDFLAGS="\$LDFLAGS -openmp" main_build_matrix_K.cpp
*
*   Microsoft Visual C++ 2008
*   mex -v COMPFLAGS="$COMPFLAGS /openmp" LINKFLAGS="$LINKFLAGS /openmp" main_build_matrix_K.cpp
*
*	Intel/{Windows}
*	mex -v CXXFLAGS="\$CXXFLAGS -Qopenmp" LDFLAGS="\$LDFLAGS -Qopenmp" main_build_matrix_K.cpp
*
*	If it doesn't support, then just use:
*	mex -v main_build_matrix_K.cpp
*/

/* How to call from Matlab:
*  [ar ai br bi] = main_build_matrix_K(nodes_glob, region_elems, region_coords, ...
*					region_node_IDs, omega(region), kappa(region), num_procs);
*  A = complex(ar,ai); B = complex(br,bi);
*  num_procs determines how many processes/cores you want to use on a single machine.
*  If num_procs > number of available CPU/cores then the value returned from
*  omp_get_num_procs() will be used.
*
*/
#include "bem.h"

#ifdef _OPENMP
#include "omp.h"
#endif

//#include "CStopWatch.h"
char first_time = 1;

#define nodes(i,j)    nodes[(i)+ nnod*(j)]
#define nod_list(i,j) nod_list[(i)+ num_nodes*(j)]
#define elements(i,j) elements[(i)+nelem*(j)]
#define elem_list(i,j) elem_list[(i)+num_elem*(j)]
#define nodes_reg(i,j) nodes_reg[(i)+nnod_reg*(j)]

#define b_real(i,j) (b_real[(i)+nnod*(j)])
#define b_img(i,j)  (b_img[(i)+nnod*(j)])

//#define _debug
/***************************************************************************
* this function is called by the mexFunction below.
* It is based on the fortran90 function of the same name.
* parameter types are declared based on that code, and according to 
* the fortran90-C conversion template provided at 
*	http://docs.hp.com/en/B3909-90002/ch08s01.html
*
* Complex numbers are however represented with the mxArray matlab data
* type, and matlab functions are used to manipulate it
*
* @param nodes: (real*8 array)
* @param elements: (integer array)
* @param nodes_reg: (real*8 2D array)
* @param inreg: (integer array)
* @param nnod: (integer)
* @param nnod_reg: (integer)
* @param nelem: (integer)
* @param omega: (complex*16)
* @param D: (real*8)
* @param mynum_procs: (integer)
* @param A: (complex*16 2D array)
* @param B: (complex*16 2D array)
*
******************************************************************************/
void main_build_matrix(double *nodes, double *elements, double *nodes_reg,
          double *inreg, long nnod, long nnod_reg, long nelem, double omegar, 
	            double omegai, double D, long mynum_procs,
	  double* Ar, double* Ai, double* Br, double* Bi) {

	/* local variables */
	int ngp = 4;
	
	long inod, j, i, el, nn;
	double *b_real, *b_img, *a_real, *a_img;
	
	double *delG = new double[nelem];
    double *normal_x = new double[nelem];
    double *normal_y = new double[nelem];
    double *normal_z = new double[nelem];
    double *eta1 = new double[ngp];
    double *eta2 = new double[ngp];
    double *w = new double[ngp];
	
	for(i=0; i<ngp; ++i)
		*(eta1+i) = *(eta2+i) = *(w+i) = 0.0;
	
	#ifdef _debug
		mexPrintf("get_gausspts_3d\n");
	#endif
	/* getting gauss points for linear triangular element */
	get_gausspts_3d(eta1, eta2, w, ngp);
	
	

#ifdef _debug
	mexPrintf("get_area_normals\n");
#endif

#if defined(_OPENMP) & defined(_debug)
	double starttime = omp_get_wtime();
#endif	
	get_area_normals(elements, nodes, nelem, nnod, delG, normal_x, normal_y, normal_z);
#if defined(_OPENMP) & defined(_debug)	
	mexPrintf("\n Area calculation took: %.8lf seconds.\n\n",omp_get_wtime()-starttime);
#endif

#ifdef _debug
	mexPrintf("mxGetPr(A)\n");
#endif
	a_real = Ar;
	b_real = Br;
	a_img =  Ai;
	b_img =  Bi;

#ifdef _OPENMP
	if (mynum_procs < omp_get_num_procs())
		omp_set_num_threads(mynum_procs);
	else
		omp_set_num_threads(omp_get_num_procs());
#ifdef _debug
	mexPrintf("\n CPUs available: %d %d\n",omp_get_num_procs(),mynum_procs);
#endif
#endif

	double xi, yi, zi;
	double xel[3], yel[3], zel[3];
	long jel[3];
	double* int_val_real;
	double* int_val_img;
	#pragma omp parallel default(none) \
	  private(nn,el,inod,j,xi,yi,zi,xel,yel,zel,jel,int_val_real, \
		  int_val_img) \
	  shared(inreg,nnod,nodes_reg,nnod_reg,nelem,elements,nodes,b_real,b_img, \
	         omegar,omegai,D,delG,eta1,eta2,w,ngp,normal_x,normal_y,normal_z,Ar,Ai,Br,Bi) 
	{

#if defined(_OPENMP) & defined(_debug)
        #pragma omp master
        mexPrintf("\n  Number of threads to be used: %d\n",omp_get_num_threads());
		#pragma omp single
		mexPrintf("  %d: thread entered the parallel regions.\n",omp_get_thread_num());
#endif
		int_val_real = new double[3];
		int_val_img  = new double[3];
		
		#pragma omp for
		for(nn = 0; nn < nnod_reg; ++nn) {
			inod = (long) inreg[nn];
			for(el = 0; el < nelem; ++el) {
				jel[0] = (long) elements(el,1);	
				jel[1] = (long) elements(el,2);
				jel[2] = (long) elements(el,3);
				xel[0] = nodes(jel[0]-1,0);
				xel[1] = nodes(jel[1]-1,0);
				xel[2] = nodes(jel[2]-1,0);
				yel[0] = nodes(jel[0]-1,1);
				yel[1] = nodes(jel[1]-1,1);
				yel[2] = nodes(jel[2]-1,1);
				zel[0] = nodes(jel[0]-1,2);
				zel[1] = nodes(jel[1]-1,2);
				zel[2] = nodes(jel[2]-1,2);
			
				if(inod == jel[0]) {
					for(j=0; j<3; *(int_val_real+j)=0.0, *(int_val_img+j)=0.0, ++j);
					singular_integrand_polar(int_val_real, int_val_img, omegar, omegai,
							D, delG[el],eta1, eta2, w, ngp);
					b_real(inod-1,jel[0]-1) += int_val_real[0];
					b_img(inod-1,jel[0]-1)  += int_val_img[0];
					b_real(inod-1,jel[1]-1) += int_val_real[1];
					b_img(inod-1,jel[1]-1)  += int_val_img[1];
					b_real(inod-1,jel[2]-1) += int_val_real[2];
					b_img(inod-1,jel[2]-1)  += int_val_img[2];
				}
				
			    else if(inod == jel[1]) {
					for(j=0; j<3; *(int_val_real+j)=0.0, *(int_val_img+j)=0.0, ++j);
					singular_integrand_polar(int_val_real, int_val_img, omegar, omegai,
							D, delG[el],eta1, eta2, w, ngp);
					b_real(inod-1,jel[0]-1) += int_val_real[2];
					b_img(inod-1,jel[0]-1)  += int_val_img[2];
					b_real(inod-1,jel[1]-1) += int_val_real[0];
					b_img(inod-1,jel[1]-1)  += int_val_img[0];
					b_real(inod-1,jel[2]-1) += int_val_real[1];
					b_img(inod-1,jel[2]-1)  += int_val_img[1];
				}
			
				else if(inod == jel[2]) {
					for(j=0; j<3; *(int_val_real+j)=0.0, *(int_val_img+j)=0.0, ++j);
					singular_integrand_polar(int_val_real, int_val_img, omegar, omegai,
							D, delG[el],eta1, eta2, w, ngp);
					b_real(inod-1,jel[0]-1) += int_val_real[1];
					b_img(inod-1,jel[0]-1)  += int_val_img[1];
					b_real(inod-1,jel[1]-1) += int_val_real[2];
					b_img(inod-1,jel[1]-1)  += int_val_img[2];
					b_real(inod-1,jel[2]-1) += int_val_real[0];
					b_img(inod-1,jel[2]-1)  += int_val_img[0];
				}
			
				else {
					xi = nodes_reg(nn,0);
					yi = nodes_reg(nn,1);
					zi = nodes_reg(nn,2); 
					compute_integral_gausspts(Ar, Ai, Br, Bi, nnod, ngp, eta1, eta2, w, inod, jel,
						xel, yel, zel, xi, yi, zi, delG[el], normal_x[el], 
						normal_y[el], normal_z[el], omegar, omegai, D);
				}	
		
			} /* end inner for loop */
		} /* end outer for loop */
		delete [] int_val_real;
		delete []  int_val_img;
    } // pragma omp parallel

	/* free temporarily used local variables */
	
	delete [] delG;
	delete [] normal_x;
	delete [] normal_y;
	delete [] normal_z;
	delete [] eta1;
	delete [] eta2;
	delete [] w;
}
static void compute_integral_gausspts(double *AAr, double *AAi, double *BBr, double *BBi,
				long size, long ngpts, 
				double *eta1, double *eta2, double *w, long ith_nod,
				long *jelem, double *xelem, double *yelem, double *zelem,
				double xx, double yy, double zz, double area, double norm_x, 
				double norm_y, double norm_z, double omegar, double omegai, double diff_coeff) {	

	double phi[4], xs, ys, zs, ri, drdn;
	long j, k, t_result;
	double z_real, z_img, gi_real, gi_img, dgdr_real, dgdr_img, dgdn_real, dgdn_img;
	double param_real, param_img, *aa_real, *aa_img, *bb_real, *bb_img;
	double temp_result, xss, yss, zss, four_PI_diff;

	param_real = omegar;
	param_img = omegai;
	aa_real = AAr + ((ith_nod-1)-size);
	aa_img  = AAi + ((ith_nod-1)-size);
	bb_real = BBr + ((ith_nod-1)-size);
	bb_img  = BBi + ((ith_nod-1)-size);
	four_PI_diff = 4*PI*diff_coeff;

	for(k = ngpts -1; k>=0; --k) {
		phi[0] = eta1[k];
		phi[1] = eta2[k];
		phi[2] = 1 - eta1[k] - eta2[k];
		
		xs = xelem[0]*phi[0] + xelem[1]*phi[1] + xelem[2]*phi[2];
		ys = yelem[0]*phi[0] + yelem[1]*phi[1] + yelem[2]*phi[2];
		zs = zelem[0]*phi[0] + zelem[1]*phi[1] + zelem[2]*phi[2];

		xss = xs-xx;
		yss = ys-yy;
		zss = zs-zz;

		ri = sqrt(xss*xss + yss*yss + zss*zss);
		drdn = (norm_x*xss + norm_y*yss + norm_z*zss)/(area*ri);

		/* z = param*ri */
		z_real = (param_real)*ri;
		z_img = (param_img)*ri;


		/* gi = exp(-z)/(4*PI*ri*diff_coeff) */
		cis((-1.0*(z_real)), (-1.0*(z_img)), gi_real, gi_img);
		(gi_real) /= four_PI_diff*ri;
		(gi_img) /= four_PI_diff*ri;

		/* dgdr = gi*((-1.0/ri) - param) */
		mult_complex(gi_real, gi_img, ((-1.0/ri)-(param_real)), 
				(-1.0*(param_img)), dgdr_real, dgdr_img);

		/* dgdn = dgdr*drdn */
		dgdn_real = (dgdr_real)*drdn;
		dgdn_img = (dgdr_img)*drdn;

		/*make use of dgdn and gi so we only have to multiply by phi[j] in
		 the inner loop...effectively hoisting the multiplications into this loop */
		temp_result = area*w[k]*diff_coeff;
		dgdn_real *= temp_result;
		dgdn_img *= temp_result;
		gi_real *= area*w[k];
		gi_img *= area*w[k];

		for(j=2; j>=0; --j) {
			t_result = jelem[j]*size;
			
			/* AA[ith_nod-1][jelem[j]-1] += phi[j]*dgdn*area*w[k]*diff_coeff */
			(*(aa_real + t_result)) +=  phi[j]*dgdn_real;
			(*(aa_img + t_result)) += phi[j]*dgdn_img;

			/* BB[ith_node-1][jelem[j]-1] += phi[j]*gi*area*w[k] */
			(*(bb_real + t_result)) += phi[j]*(gi_real);
			(*(bb_img + t_result)) += phi[j]*(gi_img);
			
		}
		 
	}
	
}
void singular_integrand_polar(double* int_val_real, double* int_val_img,
 				double omegar, double omegai, double diff_coeff,
				double area, double *gp1, double *gp2, double *wgt, long ngpts) {
   
    double basis[3], tau1, tau2;
    long k, j;
    double real_part, img_part;

    double sqrttwo = sqrt(2.0);

    for(k=0; k<ngpts; ++k) {
        tau1 = sqrttwo*gp1[k]*cos(PI*gp2[k]/4.);
        tau2 = sqrttwo*gp1[k]*sin(PI*gp2[k]/4.);
        basis[0] = tau1;
        basis[1] = tau2;
        basis[2] = 1 - tau1 - tau2;

	 	/* do int_value[j] = int_value[j] + 
				exp(-param*gp1[k]*sqrt(2.0))*basis[j]*area*wgt[k]; */
		/* in two states, do the above statements. First do the exponentiation */
		real_part = omegar;
		img_part  = omegai;
		real_part = (real_part)*-1*gp1[k]*sqrttwo;
		img_part  = (img_part )*-1*gp1[k]*sqrttwo;

			/*exp(-param*pg1[k]*sqrt(2.0))*/
		cis(real_part, img_part, real_part, img_part);	
		 
        for(j=0; j<3; ++j) {
			/* then do multiplication and addition */
			(*(int_val_real+j)) += (real_part)*basis[j]*area*wgt[k];
			(*(int_val_img+j))  += ( img_part)*basis[j]*area*wgt[k];
        }
    }
		
    for(j=0; j<3; ++j) {
		/*do int_value[j] = int_value[j]*sqrt(2.0)/(16*diff_coeff) */
		(*(int_val_real+j)) *= sqrttwo/(16.*diff_coeff);
		(*(int_val_img+j))  *= sqrttwo/(16.*diff_coeff);
    }
	
}

/******************************************************************************
 *get_gausspts_3d()
 * sets the gauss point values of the 3D triangular element by modifying the 
 * formal parameter values passed in as pointers
 * derived from the f90 function of the same name
 *
 * @param eta1: (real*8 array)
 * @param eta2: (real*8 array)
 * @param w: (real*8 array)
 * @param ngp: (integer)
 *
 *****************************************************************************/
void get_gausspts_3d(double *eta1, double *eta2, double *w, int ngp) {
	
	/* getting gass points for 3D triangular element */
	if(ngp == 3) {
		eta1[0] = 1.0/2.0;
		eta1[1] = 0.0;
		eta1[2] = 1.0/2.0;
		eta2[0] = 1.0/2.0;
		eta2[1] = 1.0/2.0;
		eta2[2] = 0.0;
		
		w[0]    = 1.0/3.0;
		w[1]    = 1.0/3.0;
		w[2]    = 1.0/3.0;

	}
	else 
		if(ngp == 4) {
			eta1[0] = 1.0/3.0;
			eta1[1] = 3.0/5.0;
			eta1[2] = 1.0/5.0;
			eta1[3] = 1.0/5.0;
			eta2[0] = 1.0/3.0;
			eta2[1] = 1.0/5.0;
			eta2[2] = 3.0/5.0;
			eta2[3] = 1.0/5.0;
		
			w[0]    = -27.0/48.0;
			w[1]    =  25.0/48.0;
			w[2]    =  25.0/48.0;
			w[3]    =  25.0/48.0;
		}
}

/*************************************************************************************
 * get_area_normals()
 * based on the f90 subroutine of the same name
 *
 * @param elem_list: (integer 2D array)
 * @param nod_list: (real*8 2D array)
 * @param num_elem: (integer)
 * @param num_nodes: (integer)
 * @param area_elem_list: (real*8 2D array)
 * @param normal_x_list: (real*8 array)
 * @param normal_y_list: (real*8 array)
 * @param normal_z_list: (real*8 array)
 *
 ************************************************************************************/
void get_area_normals(double *elem_list, double *nod_list, long num_elem,
			long num_nodes, double *area_elem_list, double *normal_x_list,
			double *normal_y_list, double *normal_z_list) {

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double dx_deta1, dx_deta2, dy_deta1, dy_deta2, dz_deta1, dz_deta2;
	long jelem[3], ee;

	#pragma omp parallel default(none) \
		shared(normal_x_list,normal_y_list,normal_z_list,area_elem_list, \
			   nod_list,elem_list,num_elem,num_nodes) \
		private(ee,jelem,x1,x2,x3,y1,y2,y3,z1,z2,z3,dx_deta1,dx_deta2, \
		        dy_deta1,dy_deta2,dz_deta1,dz_deta2)
	{
		#pragma omp for
		for(ee = 0; ee < num_elem; ++ee) {
			jelem[0] = (long) elem_list(ee,1);
			jelem[1] = (long) elem_list(ee,2);
			jelem[2] = (long) elem_list(ee,3);

			x1 = nod_list(jelem[0]-1,0);
			x2 = nod_list(jelem[1]-1,0);
			x3 = nod_list(jelem[2]-1,0);

			y1 = nod_list(jelem[0]-1,1);
			y2 = nod_list(jelem[1]-1,1);
			y3 = nod_list(jelem[2]-1,1);

			z1 = nod_list(jelem[0]-1,2);
			z2 = nod_list(jelem[1]-1,2);
			z3 = nod_list(jelem[2]-1,2);

			dx_deta1 = x1 - x3;
			dx_deta2 = x2 - x3;
			dy_deta1 = y1 - y3;
			dy_deta2 = y2 - y3;
			dz_deta1 = z1 - z3;
			dz_deta2 = z2 - z3; 

			normal_x_list[ee] = dy_deta1*dz_deta2 - dy_deta2*dz_deta1;
			normal_y_list[ee] = dz_deta1*dx_deta2 - dx_deta1*dz_deta2;
			normal_z_list[ee] = dx_deta1*dy_deta2 - dy_deta1*dx_deta2; 
			area_elem_list[ee] = sqrt(pow(normal_x_list[ee],2) +
									pow(normal_y_list[ee],2) + pow(normal_z_list[ee],2));
		}
	}
}
/***************************************************************************************
* mexFunction()
* gateway from matlab workspace.
* based on the f90 mexFunction for main_build_matrix.f90
*
* @param nlhs: number of left hand side return value (returned by reference)
* @param plhs: array of pointers to l.h.s formal parameters
* @param nrhs: number of right hand side parameters
* @param prhs: array of pointers to r.h.s formal parameter
*
****************************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	long nnod, nelem, nnod_reg, mynum_procs;
	double *elements;
	double *inreg;
	double *nodes, *nodes_reg; /* simulate 2D array */
	double D;
	double omegar, omegai, *AAr, *AAi, *BBr, *BBi;

	
	/* parameter checking. ensure 7 inputs and 4 outputs vars */
	if(nrhs != 7) {
		mexErrMsgTxt("main_build_matrix_K.mex: 7 inputs are required\n");
	} 
	if(nlhs != 4)
		mexErrMsgTxt("main_build_matrix_K.mex: 4 outputs required: [Ar Ai Br Bi] = main_matrix...\n");

	/* nodes matrix */
	nnod = mxGetM(prhs[0]);
	nodes = mxGetPr(prhs[0]);
	
	/* element matrix */
	nelem = mxGetM(prhs[1]);
	elements = mxGetPr(prhs[1]);

	/* nodereg matrix */
	nnod_reg = mxGetM(prhs[2]);
	nodes_reg = mxGetPr(prhs[2]);

	/* inreg array */
	inreg = (double *) mxGetData(prhs[3]);
	
	/* copy omega */
	omegar = *(mxGetPr(prhs[4]));
	if (mxIsComplex(prhs[4]))
		omegai = *(mxGetPi(prhs[4]));
	else {
		omegai = 0.;
	}
	
	/* D */
	D = *(mxGetPr(prhs[5]));


	/* Number of cores to be used */
	mynum_procs = (long) *(mxGetPr(prhs[6]));
	
	plhs[0] = mxCreateDoubleMatrix(nnod,nnod,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(nnod,nnod,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(nnod,nnod,mxREAL);
	plhs[3] = mxCreateDoubleMatrix(nnod,nnod,mxREAL);
	
	AAr = mxGetPr(plhs[0]);
	AAi = mxGetPr(plhs[1]);
	BBr = mxGetPr(plhs[2]);
	BBi = mxGetPr(plhs[3]);

	/* call the computational routines */
#if defined(_OPENMP) & defined(_debug)
	double starttime = omp_get_wtime();
#endif

	main_build_matrix(nodes, elements, nodes_reg, inreg, nnod, nnod_reg, nelem, \
		              omegar, omegai, D, mynum_procs, AAr, AAi, BBr, BBi);
		
#if defined(_OPENMP) & defined(_debug)
	double endtime = omp_get_wtime();
	mexPrintf(" Elapsed time: %.8lf seconds.\n\n",endtime-starttime);
#endif
}
