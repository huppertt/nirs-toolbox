#include <stdio.h>
#include <math.h>
#include <mex.h>

/* Integrates light source at each node
 *
 * Part of NIRFAST package
 * H Dehghani
 * Last Updated - 7/13/09 M Jermyn
 */

double phidotphi_tri3(double g[3][2], double mua[3], double val[3][3]);
double phidotphi_tet4(double g[4][3], double mua[4], double val[4][4]);

/* -------- Heart of the mex file----------- */
void mainloop(double *nodes,
	      double *elements,
	      double *dim,
	      double *source,
	      double *width,
	      int nodem,
	      int noden,
	      int elemm,
	      int elemn,
	      double *Sr,
	      double *Si)
{
  double w2;
  #define pi 3.14159265358979

  w2 = *width;

  if (dim[0] == 2){
    int ele, i, j, index[3], ind;
    double indextmp, d, q, tri_vtx[3][2];
    double disttmp[3], c_val[3][3], k[3];
    static double ones[3] = {1, 1, 1};
    double *dist;
    dist=(double *)mxMalloc((size_t) (nodem*sizeof(double)));

    for (i=0; i<nodem; ++i){
      dist[i] = sqrt((*(nodes+i) - source[0]) *	      \
		     (*(nodes+i) - source[0]) +	      \
		     (*(nodes+(i+nodem)) - source[1]) *     \
		     (*(nodes+(i+nodem)) - source[1]));
    }

    for (ele=0; ele<elemm; ++ele){
      for (i=0; i<elemn; ++i){
	indextmp = *(elements+(ele+(i*elemm)));
	index[i] = indextmp;
	k[i] = indextmp;
	disttmp[i] = *(dist+(index[i]-1));
      }
      if (disttmp[0] <= w2 && disttmp[1] <= w2 && disttmp[2] <= w2){
	for (i=0; i<elemn; ++i){
	  for (j=0; j<noden; ++j){
	    tri_vtx[i][j] = *(nodes+(index[i]-1+(j*nodem)));
	  }
	}
      	for (i=0; i<3; ++i){
	  d = disttmp[i];
	  q = exp (-0.5 * d*d) / sqrt(2.0*pi);
	  phidotphi_tri3(tri_vtx,ones,c_val);
	  for (j=0; j<3; ++j){
	    ind = *(k+j);
	    Sr[ind-1] += q*c_val[i][j]*1*cos(0.15);
	    Si[ind-1] += q*c_val[i][j]*1*sin(0.15);
	  }
	}	  
      }
    }
  }

  else if (dim[0] == 3){
    int ele, i, j, index[4], ind;
    double indextmp, d, q, tri_vtx[4][3];
    double disttmp[4], c_val[4][4], k[4];
    static double ones[4] = {1, 1, 1, 1};
    double *dist;
    dist=(double *)mxMalloc((size_t) (nodem*sizeof(double)));

    for (i=0; i<nodem; ++i){
      dist[i] = sqrt((*(nodes+i) - source[0]) *	      \
		     (*(nodes+i) - source[0]) +	      \
		     (*(nodes+(i+nodem)) - source[1]) *     \
		     (*(nodes+(i+nodem)) - source[1]) +	    \
		     (*(nodes+(i+nodem+nodem)) - source[2]) *	\
		     (*(nodes+(i+nodem+nodem)) - source[2]));
    }
    for (ele=0; ele<elemm; ++ele){
      for (i=0; i<elemn; ++i){
	indextmp = *(elements+(ele+(i*elemm)));
	index[i] = indextmp;
	k[i] = indextmp;
	disttmp[i] = *(dist+(index[i]-1));
      }
      if (disttmp[0] <= w2 && disttmp[1] <= w2 && \
	  disttmp[2] <= w2 && disttmp[3] <= w2){
	for (i=0; i<elemn; ++i){
	  for (j=0; j<noden; ++j){
	    tri_vtx[i][j] = *(nodes+(index[i]-1+(j*nodem)));
	  }
	}
      	for (i=0; i<4; ++i){
	  d = disttmp[i];
	  q = exp (-0.5 * d*d) / sqrt(2.0*pi);
	  phidotphi_tet4(tri_vtx,ones,c_val);
	  for (j=0; j<4; ++j){
	    ind = *(k+j);
	    Sr[ind-1] += q*c_val[i][j]*1*cos(0.15);
	    Si[ind-1] += q*c_val[i][j]*1*sin(0.15);
	  }
	}	  
      }
    }
  }

  return;
}


double phidotphi_tri3(double g[3][2], double mua[3], double val[3][3])
{
  int ii, jj;
  static int L[2][3] = {{-1.0,1.0,0.0}, \
			{-1.0,0.0,1.0}};
  static double w[3] = {1.0/6.0,1.0/6.0,1.0/6.0};
  static double ip[3][2] = {{1.0/2.0,0.0},{1.0/2.0,1.0/2.0},{0.0,1.0/2.0}};
  double Jt[2][2], dJt, S[3], tmp;
  double x0 = g[0][0]; double x1 = g[1][0];
  double x2 = g[2][0];
  double y0 = g[0][1]; double y1 = g[1][1];
  double y2 = g[2][1];
  double k0 = mua[0]; double k1 = mua[1];
  double k2 = mua[2];

  Jt[0][0] = L[0][0]*x0 + L[0][1]*x1 + L[0][2]*x2;
  Jt[0][1] = L[0][0]*y0 + L[0][1]*y1 + L[0][2]*y2;
  Jt[1][0] = L[1][0]*x0 + L[1][1]*x1 + L[1][2]*x2;
  Jt[1][1] = L[1][0]*y0 + L[1][1]*y1 + L[1][2]*y2;
  
  dJt = (Jt[0][0]*Jt[1][1]) - (Jt[0][1]*Jt[1][0]);
  
  dJt = sqrt(dJt*dJt);
    
  for (ii=0; ii<3; ++ii){
    for (jj=0; jj<3; ++jj){
      val[ii][jj]=0;
    }
  }

  for (ii=0; ii<3; ++ii){
    S[0] = 1-ip[ii][0]-ip[ii][1];
    S[1] = ip[ii][0]; 
    S[2] = ip[ii][1];
    tmp = k0*S[0] + k1*S[1] + k2*S[2];
    tmp = tmp*w[ii];
    val[0][0] += S[0]*S[0]*tmp;
    val[0][1] += S[1]*S[0]*tmp;
    val[0][2] += S[2]*S[0]*tmp;
    val[1][0] += S[0]*S[1]*tmp;
    val[1][1] += S[1]*S[1]*tmp;
    val[1][2] += S[2]*S[1]*tmp;
    val[2][0] += S[0]*S[2]*tmp;
    val[2][1] += S[1]*S[2]*tmp;
    val[2][2] += S[2]*S[2]*tmp;
    }
  for (ii=0; ii<3; ++ii){
    for (jj=0; jj<3; ++jj){
      val[ii][jj] = val[ii][jj]*dJt;
    }
  }
  return;
}

double phidotphi_tet4(double g[4][3], double mua[4], double val[4][4])
{
  int ii, jj;
  static int L[3][4] = {{-1.0,1.0,0.0,0.0}, \
			{-1.0,0.0,1.0,0.0}, \
			{-1.0,0.0,0.0,1.0}};
  static double w[24] = {6.6537917096945820e-3, 6.6537917096945820e-3, \
			 6.6537917096945820e-3, 6.6537917096945820e-3, \
			 1.6795351758867738e-3, 1.6795351758867738e-3, \
			 1.6795351758867738e-3, 1.6795351758867738e-3, \
			 9.2261969239424536e-3, 9.2261969239424536e-3, \
			 9.2261969239424536e-3, 9.2261969239424536e-3, \
			 8.0357142857142857e-3, 8.0357142857142857e-3, \
			 8.0357142857142857e-3, 8.0357142857142857e-3, \
			 8.0357142857142857e-3, 8.0357142857142857e-3, \
			 8.0357142857142857e-3, 8.0357142857142857e-3, \
			 8.0357142857142857e-3, 8.0357142857142857e-3, \
			 8.0357142857142857e-3, 8.0357142857142857e-3};
  static double ip[24][3] = \
    {{0.21460287125915202,0.21460287125915202,0.3561913862225439}, \
     {0.21460287125915202,0.3561913862225439,0.21460287125915202}, \
     {0.3561913862225439,0.21460287125915202,0.21460287125915202}, \
     {0.21460287125915202,0.21460287125915202,0.21460287125915202}, \
     {0.040673958534611353,0.040673958534611353,0.87797812439616594}, \
     {0.040673958534611353,0.87797812439616594,0.040673958534611353}, \
     {0.87797812439616594,0.040673958534611353,0.040673958534611353}, \
     {0.040673958534611353,0.040673958534611353,0.040673958534611353}, \
     {0.32233789014227551,0.32233789014227551,0.0329863295731735}, \
     {0.32233789014227551,0.0329863295731735,0.32233789014227551}, \
     {0.0329863295731735,0.32233789014227551,0.32233789014227551}, \
     {0.32233789014227551,0.32233789014227551,0.32233789014227551}, \
     {0.063661001875017525,0.26967233145831580,0.6030056647916492}, \
     {0.063661001875017525,0.6030056647916492,0.26967233145831580}, \
     {0.26967233145831580,0.063661001875017525,0.6030056647916492}, \
     {0.26967233145831580,0.6030056647916492,0.063661001875017525}, \
     {0.6030056647916492,0.063661001875017525,0.26967233145831580}, \
     {0.6030056647916492,0.26967233145831580,0.063661001875017525}, \
     {0.063661001875017525,0.063661001875017525,0.6030056647916492}, \
     {0.063661001875017525,0.6030056647916492,0.063661001875017525}, \
     {0.6030056647916492,0.063661001875017525,0.063661001875017525}, \
     {0.063661001875017525,0.063661001875017525,0.26967233145831580}, \
     {0.063661001875017525,0.26967233145831580,0.063661001875017525}, \
     {0.26967233145831580,0.063661001875017525,0.063661001875017525}};
  double Jt[3][3], dJt, S[4][1], tmp;
  double x0 = g[0][0]; double x1 = g[1][0];
  double x2 = g[2][0]; double x3 = g[3][0];
  double y0 = g[0][1]; double y1 = g[1][1];
  double y2 = g[2][1]; double y3 = g[3][1];
  double z0 = g[0][2]; double z1 = g[1][2];
  double z2 = g[2][2]; double z3 = g[3][2];
  double k0 = mua[0]; double k1 = mua[1];
  double k2 = mua[2]; double k3 = mua[3];
  
  Jt[0][0] = L[0][0]*x0 + L[0][1]*x1 + L[0][2]*x2 + L[0][3]*x3;
  Jt[0][1] = L[0][0]*y0 + L[0][1]*y1 + L[0][2]*y2 + L[0][3]*y3;
  Jt[0][2] = L[0][0]*z0 + L[0][1]*z1 + L[0][2]*z2 + L[0][3]*z3;
  Jt[1][0] = L[1][0]*x0 + L[1][1]*x1 + L[1][2]*x2 + L[1][3]*x3;
  Jt[1][1] = L[1][0]*y0 + L[1][1]*y1 + L[1][2]*y2 + L[1][3]*y3;
  Jt[1][2] = L[1][0]*z0 + L[1][1]*z1 + L[1][2]*z2 + L[1][3]*z3;
  Jt[2][0] = L[2][0]*x0 + L[2][1]*x1 + L[2][2]*x2 + L[2][3]*x3;
  Jt[2][1] = L[2][0]*y0 + L[2][1]*y1 + L[2][2]*y2 + L[2][3]*y3;
  Jt[2][2] = L[2][0]*z0 + L[2][1]*z1 + L[2][2]*z2 + L[2][3]*z3;
  
  dJt = Jt[0][0]*(Jt[1][1]*Jt[2][2] - Jt[1][2]*Jt[2][1]) \
    - Jt[1][0]*(Jt[0][1]*Jt[2][2] - Jt[0][2]*Jt[2][1]) \
    + Jt[2][0]*(Jt[0][1]*Jt[1][2] - Jt[0][2]*Jt[1][1]);
  
  dJt = sqrt(dJt*dJt);
  
  for (ii=0; ii<4; ++ii){
    for (jj=0; jj<4; ++jj){
      val[ii][jj]=0;
    }
  }

  for (ii=0; ii<24; ++ii){
    S[0][0] = 1-ip[ii][0]-ip[ii][1]-ip[ii][2];
    S[1][0] = ip[ii][0]; 
    S[2][0] = ip[ii][1];
    S[3][0] = ip[ii][2];
    tmp = k0*S[0][0] + k1*S[1][0] + k2*S[2][0] + k3*S[3][0];
    tmp = tmp*w[ii];
    val[0][0] += S[0][0]*S[0][0]*tmp;
    val[0][1] += S[1][0]*S[0][0]*tmp;
    val[0][2] += S[2][0]*S[0][0]*tmp;
    val[0][3] += S[3][0]*S[0][0]*tmp;
    val[1][0] += S[0][0]*S[1][0]*tmp;
    val[1][1] += S[1][0]*S[1][0]*tmp;
    val[1][2] += S[2][0]*S[1][0]*tmp;
    val[1][3] += S[3][0]*S[1][0]*tmp;
    val[2][0] += S[0][0]*S[2][0]*tmp;
    val[2][1] += S[1][0]*S[2][0]*tmp;
    val[2][2] += S[2][0]*S[2][0]*tmp;
    val[2][3] += S[3][0]*S[2][0]*tmp;
    val[3][0] += S[0][0]*S[3][0]*tmp;
    val[3][1] += S[1][0]*S[3][0]*tmp;
    val[3][2] += S[2][0]*S[3][0]*tmp;
    val[3][3] += S[3][0]*S[3][0]*tmp;
    }
  for (ii=0; ii<4; ++ii){
    for (jj=0; jj<4; ++jj){
      val[ii][jj] = val[ii][jj]*dJt;
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
  double *nodes,*elements,*source, *width, *dim;
  int nodem,noden,elemm,elemn;
  double *Sr,*Si;
  
  /* Error checking  */
  
  if (nrhs < 5 )
    mexErrMsgTxt(" There is not enough input arguments");
  
  if(nlhs!=1)
    mexErrMsgTxt("This routine requires One ouput argument");
  
  nodes=mxGetPr(prhs[0]);          /*  nodes of mesh */
  elements=mxGetPr(prhs[1]);       /*  elements of mesh */
  dim=mxGetPr(prhs[2]);            /*  mesh.dimension */
  source=mxGetPr(prhs[3]);         /*  source location */
  width=mxGetPr(prhs[4]);          /*  width of source */

  nodem=mxGetM(prhs[0]);          /*  Number of of nodes */
  noden=mxGetN(prhs[0]);          /*  Number of node freedom */
  elemm=mxGetM(prhs[1]);          /*  Number of elements */
  elemn=mxGetN(prhs[1]);          /*  Number of nodes per elements */

  plhs[0]=mxCreateDoubleMatrix(nodem,1,mxCOMPLEX); /* s */
  
  Sr=mxGetPr(plhs[0]);
  Si=mxGetPi(plhs[0]);
  
  mainloop(nodes,elements,dim,source,width,nodem,noden,elemm,elemn,Sr,Si);
  
  return;
}
