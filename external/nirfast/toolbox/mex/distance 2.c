#include <stdio.h>
#include <math.h>
#include <mex.h>

/* Calculates the distance of a node from every boundary node
 *
 * Part of NIRFAST package
 * H Dehghani
 * Last Updated - 7/13/09 M Jermyn
 */

/* -------- Heart of the mex file----------- */
void mainloop(double *nodes,
	      double *bndvtx,
	      double *point,
	      double *dist,
	      int nodem)
{
  int i; 

  for (i=0; i<nodem; ++i){
    if (bndvtx[i] == 1){
      dist[i] = sqrt((*(nodes+i) - point[0]) * (*(nodes+i) - point[0]) + (*(nodes+(i+nodem)) - point[1]) * (*(nodes+(i+nodem)) - point[1]) + (*(nodes+(i+nodem+nodem)) - point[2]) * (*(nodes+(i+nodem+nodem)) - point[2]));
    }
    else {
      dist[i] = 1000;
    }
  }
  return;
}
/* -------- Gate-way to matlab  ------------ */

void mexFunction(int nlhs,
		 mxArray *plhs[],
		 int nrhs,
		 const mxArray *prhs[])

{
  double *nodes,*bndvtx,*point;
  int nodem;
  double *dist;
  
  /* Error checking  */

  if (nrhs < 2 )
    mexErrMsgTxt(" There is not enough input arguments");
  
  if(nlhs!=1)
    mexErrMsgTxt("This routine requires one ouput arguments");
  
  nodes=mxGetPr(prhs[0]);          /*  nodes of mesh */
  bndvtx=mxGetPr(prhs[1]);         /*  boundary nodes of mesh */
  point=mxGetPr(prhs[2]);          /*  boundary point */
  
  nodem=mxGetM(prhs[0]);          /*  Number of rows of nodes */
  
  plhs[0]=mxCreateDoubleMatrix(nodem,1,mxREAL); /* vector to return to Matlab */
  
  dist=mxGetPr(plhs[0]);
  
  mainloop(nodes,bndvtx,point,dist,nodem);
  
  return;
}
