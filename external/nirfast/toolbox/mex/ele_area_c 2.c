#include <stdio.h>
#include <math.h>
#include <mex.h>

/* Calculates area of elements
 *
 * Part of NIRFAST package
 * H Dehghani
 * Last Updated - 7/13/09 M Jermyn
 */

/* -------- Heart of the mex file----------- */
void mainloop(double *nodes,
	      double *elements,
	      int nodem,
	      int noden,
	      int elemm,
	      int elemn,
	      double *I)
{
  if (elemn == 3){
    int i, j, k, ele, index[3];
    double indextmp, tri_vtx[3][2];
    static int L[2][3] = {{-1.0,1.0,0.0},{-1.0,0.0,1.0}};
    double Jt[2][2], dJt;
    k = 0;
    for (ele=0; ele<elemm; ++ele){
      for (i=0; i<elemn; ++i){
	indextmp = *(elements+(ele+(i*elemm)));
	index[i] = indextmp;
      }
      for (i=0; i<elemn; ++i){
	for (j=0; j<noden; ++j){
	  tri_vtx[i][j] = *(nodes+(index[i]-1+(j*nodem)));
	}
      }

      Jt[0][0] = L[0][0]*tri_vtx[0][0] + L[0][1]*tri_vtx[1][0] + L[0][2]*tri_vtx[2][0];
      Jt[0][1] = L[0][0]*tri_vtx[0][1] + L[0][1]*tri_vtx[1][1] + L[0][2]*tri_vtx[2][1];
      Jt[1][0] = L[1][0]*tri_vtx[0][0] + L[1][1]*tri_vtx[1][0] + L[1][2]*tri_vtx[2][0];
      Jt[1][1] = L[1][0]*tri_vtx[0][1] + L[1][1]*tri_vtx[1][1] + L[1][2]*tri_vtx[2][1];
      
      dJt = Jt[0][0]*Jt[1][1] - Jt[0][1]*Jt[1][0];
      
      dJt = sqrt(dJt*dJt);
      
      I[k] = dJt*0.5;
      ++k;
    }    
  }
  else if (elemn == 4){
    int i, j, k, ele, index[4];
    double indextmp, tri_vtx[4][3];
    static int L[3][4] = {{-1.0,1.0,0.0,0.0},{-1.0,0.0,1.0,0.0},{-1.0,0.0,0.0,1.0}};
    double Jt[3][3], dJt;
    k = 0;
    for (ele=0; ele<elemm; ++ele){
      for (i=0; i<elemn; ++i){
	indextmp = *(elements+(ele+(i*elemm)));
	index[i] = indextmp;
      }
      for (i=0; i<elemn; ++i){
	for (j=0; j<noden; ++j){
	  tri_vtx[i][j] = *(nodes+(index[i]-1+(j*nodem)));
	}
      }

      Jt[0][0] = L[0][0]*tri_vtx[0][0] + L[0][1]*tri_vtx[1][0] + L[0][2]*tri_vtx[2][0] + L[0][3]*tri_vtx[3][0];
      Jt[0][1] = L[0][0]*tri_vtx[0][1] + L[0][1]*tri_vtx[1][1] + L[0][2]*tri_vtx[2][1] + L[0][3]*tri_vtx[3][1];
      Jt[0][2] = L[0][0]*tri_vtx[0][2] + L[0][1]*tri_vtx[1][2] + L[0][2]*tri_vtx[2][2] + L[0][3]*tri_vtx[3][2];
      Jt[1][0] = L[1][0]*tri_vtx[0][0] + L[1][1]*tri_vtx[1][0] + L[1][2]*tri_vtx[2][0] + L[1][3]*tri_vtx[3][0];
      Jt[1][1] = L[1][0]*tri_vtx[0][1] + L[1][1]*tri_vtx[1][1] + L[1][2]*tri_vtx[2][1] + L[1][3]*tri_vtx[3][1];
      Jt[1][2] = L[1][0]*tri_vtx[0][2] + L[1][1]*tri_vtx[1][2] + L[1][2]*tri_vtx[2][2] + L[1][3]*tri_vtx[3][2];
      Jt[2][0] = L[2][0]*tri_vtx[0][0] + L[2][1]*tri_vtx[1][0] + L[2][2]*tri_vtx[2][0] + L[2][3]*tri_vtx[3][0];
      Jt[2][1] = L[2][0]*tri_vtx[0][1] + L[2][1]*tri_vtx[1][1] + L[2][2]*tri_vtx[2][1] + L[2][3]*tri_vtx[3][1];
      Jt[2][2] = L[2][0]*tri_vtx[0][2] + L[2][1]*tri_vtx[1][2] + L[2][2]*tri_vtx[2][2] + L[2][3]*tri_vtx[3][2];
      
      dJt = Jt[0][0]*(Jt[1][1]*Jt[2][2] - Jt[1][2]*Jt[2][1]) \
        - Jt[1][0]*(Jt[0][1]*Jt[2][2] - Jt[0][2]*Jt[2][1]) \
        + Jt[2][0]*(Jt[0][1]*Jt[1][2] - Jt[0][2]*Jt[1][1]);
      
      dJt = sqrt(dJt*dJt);

      I[k] = dJt/6;
      ++k;
    }
  }
}
/* -------- Gate-way to matlab  ------------ */
void mexFunction(int nlhs,
		 mxArray *plhs[],
		 int nrhs,
		 const mxArray *prhs[])
     
{
  double *nodes,*elements;
  int nodem,noden,elemm,elemn;
  double *I;
  
  /* Error checking  */
  
  if (nrhs < 2 )
    mexErrMsgTxt(" There is not enough input arguments");
  
  if(nlhs!=1)
    mexErrMsgTxt("This routine requires 1 ouput arguments");
  
  nodes=mxGetPr(prhs[0]);          /*  nodes of mesh */
  elements=mxGetPr(prhs[1]);       /*  elements of mesh */
  
  nodem=mxGetM(prhs[0]);          /*  Number of of nodes */
  noden=mxGetN(prhs[0]);          /*  Number of node freedom */
  elemm=mxGetM(prhs[1]);          /*  Number of elements */
  elemn=mxGetN(prhs[1]);          /*  Number of nodes per elements */
  
  plhs[0]=mxCreateDoubleMatrix(elemm,1,mxREAL);    /* I */

  I=mxGetPr(plhs[0]);
  
  mainloop(nodes,elements,nodem,noden,elemm,elemn,I);
  
  return;
}
