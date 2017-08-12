 #include <stdio.h>
#include <math.h>
#include <mex.h>

/* Sums all element areas around a node
 *
 * Part of NIRFAST package
 * H Dehghani
 * Last Updated - 7/13/09 M Jermyn
 */

/* -------- Heart of the mex file----------- */
void mainloop(double *nodes,
	      double *elements,
	      double *area,
	      int nodem,
	      int noden,
	      int elemm,
	      int elemn,
	      double *support)

{
  int i, j, k;

  for (i=0; i<elemm; ++i){
    for (j=0; j<(elemn); ++j){
      k = *(elements+(i+(j*elemm)));
      support[k-1]  = support[k-1] + *(area+(i));
    }
  }
}

/* -------- Gate-way to matlab  ------------ */
void mexFunction(int nlhs,
		 mxArray *plhs[],
		 int nrhs,
		 const mxArray *prhs[])
     
{
  double *nodes,*elements,*area;
  int nodem,noden,elemm,elemn;
  double *support;
  
  /* Error checking  */
  
  if (nrhs < 3 )
    mexErrMsgTxt(" There is not enough input arguments");
  
  if(nlhs!=1)
    mexErrMsgTxt("This routine requires 1 ouput arguments");
  
  nodes=mxGetPr(prhs[0]);          /*  nodes of mesh */
  elements=mxGetPr(prhs[1]);       /*  elements of mesh */
  area=mxGetPr(prhs[2]);           /*  element area */

  nodem=mxGetM(prhs[0]);          /*  Number of of nodes */
  noden=mxGetN(prhs[0]);          /*  Number of node freedom */
  elemm=mxGetM(prhs[1]);          /*  Number of elements */
  elemn=mxGetN(prhs[1]);          /*  Number of nodes per elements */
  
  plhs[0]=mxCreateDoubleMatrix(nodem,1,mxREAL);    /* support */

  support=mxGetPr(plhs[0]);
  
  mainloop(nodes,elements,area,nodem,noden,elemm,elemn,support);
  
  return;
}
