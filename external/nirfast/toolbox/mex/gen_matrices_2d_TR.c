#include <stdio.h>
#include <math.h>
#include <mex.h>
double gradphi(double g[3][2], double kappa[3], double val[3][3]);
double phidotphi(double g[3][2], double mua[3], double val[3][3]);
double boundary_int(double g[3][2], double bnd_index[3],	\
		    double ksir[3], double valr[2][2]);
double bound2(double g[2][2], int il, int im, double ksir[2], double *val);

/* -------- Heart of the mex file----------- */
void mainloop(double *nodes,
	      double *elements,
	      double *bndvtx,
	      double *mua,
	      double *kappa,
	      double *ksir,
	      int nodem,
	      int noden,
	      int elemm,
	      int elemn,
          double *dt,
	      double *I,
	      double *J,
	      double *S1,
	      double *S2)
{
  int ele, i, j, index[3], k;
  double indextmp, tri_vtx[3][2];
  double bnd_index[3], kappa_index[3], mua_index[3], dt_index[3];
  double k_val[3][3], c_val[3][3], f_valr[2][2], m_val[3][3];
  double ksir_index[3], unity[3], dtrec[1];
  
  for (i=0; i<3; ++i){
    unity[i]= 1.0;
  }
  
  /*dtrec[0] = 1/(*dt);*/
  /*mexPrintf("%g\n",dtrec[0]);*/
  
  k = 0;
  for (ele=0; ele<elemm; ++ele){
    for (i=0; i<elemn; ++i){
      indextmp = *(elements+(ele+(i*elemm)));
      index[i] = indextmp;
    }

    for (i=0; i<elemn; ++i){
      kappa_index[i] = *(kappa+(index[i]-1));
      mua_index[i] = *(mua+(index[i]-1));
      bnd_index[i] = *(bndvtx+(index[i]-1));
      ksir_index[i] = *(ksir+(index[i]-1));
      dt_index[i] = 1/(*(dt+(index[i]-1)));
     for (j=0; j<noden; ++j){
	tri_vtx[i][j] = *(nodes+(index[i]-1+(j*nodem)));
      }
    }
    gradphi(tri_vtx,kappa_index,k_val);
    phidotphi(tri_vtx,mua_index,c_val);
    phidotphi(tri_vtx,unity,m_val);

    for (i=0; i<3; ++i){
     for (j=0; j<3; ++j){
       I[k] = index[i];
       J[k] = index[j];
       S1[k] = (0.5*k_val[j][i])+(0.5*c_val[j][i])+(dt_index[i]*m_val[j][i]);
       S2[k] = (0.5*k_val[j][i])+(0.5*c_val[j][i])-(dt_index[i]*m_val[j][i]);
      ++k;
     }
    }


    if (bnd_index[0]+bnd_index[1]+bnd_index[2] != 0){
      boundary_int(tri_vtx, bnd_index, ksir_index, f_valr);
      if (bnd_index[0]==1 && bnd_index[2]==1){
	int index1[2];
	index1[0] = index[0];
	index1[1] = index[2];
	for (i=0; i<2; ++i){
	  for (j=0; j<=i; ++j){
	    I[k] = index1[i];
	    J[k] = index1[j];
	    S1[k] = 0.5*f_valr[i][j];
	    S2[k] = 0.5*f_valr[i][j];
	    ++k;
	  }
	}
      }
      else if (bnd_index[0]==1 && bnd_index[1]==1){
	int index1[2];
	index1[0] = index[0];
	index1[1] = index[1];
	for (i=0; i<2; ++i){
	  for (j=0; j<=i; ++j){
	    I[k] = index1[i];
	    J[k] = index1[j];
	    S1[k] = 0.5*f_valr[i][j];
	    S2[k] = 0.5*f_valr[i][j];
	    ++k;
	  }
	}
      }
      else if (bnd_index[1]==1 && bnd_index[2]==1){
	int index1[2];
	index1[0] = index[1];
	index1[1] = index[2];
	for (i=0; i<2; ++i){
	  for (j=0; j<=i; ++j){
	    I[k] = index1[i];
	    J[k] = index1[j];
	    S1[k] = 0.5*f_valr[i][j];
	    S2[k] = 0.5*f_valr[i][j];
	    ++k;
	  }
	}
      }
    }

 
  }
  return;
}
double gradphi(double g[3][2], double kappa[3], double val[3][3])
{
  int ii, jj;
  static int L[2][3] = {{-1.0,1.0,0.0}, \
			{-1.0,0.0,1.0}};
  static double w[3] = {1.0/6.0,1.0/6.0,1.0/6.0};
  static double ip[3][2] = {{1.0/2.0,0.0},{1.0/2.0,1.0/2.0},{0.0,1.0/2.0}};
  double Jt[2][2], dJt, iJt[2][2], S[3], G[2][3], GG[3][3], tmp;
  double x0 = g[0][0]; double x1 = g[1][0];
  double x2 = g[2][0];
  double y0 = g[0][1]; double y1 = g[1][1];
  double y2 = g[2][1];
  double k0 = kappa[0]; double k1 = kappa[1];
  double k2 = kappa[2];
  
  Jt[0][0] = L[0][0]*x0 + L[0][1]*x1 + L[0][2]*x2;
  Jt[0][1] = L[0][0]*y0 + L[0][1]*y1 + L[0][2]*y2;
  Jt[1][0] = L[1][0]*x0 + L[1][1]*x1 + L[1][2]*x2;
  Jt[1][1] = L[1][0]*y0 + L[1][1]*y1 + L[1][2]*y2;
  
  dJt = (Jt[0][0]*Jt[1][1]) - (Jt[0][1]*Jt[1][0]);
  
  iJt[0][0] = +(Jt[1][1])/dJt;
  iJt[0][1] = -(Jt[0][1])/dJt;
  iJt[1][0] = -(Jt[1][0])/dJt;
  iJt[1][1] = +(Jt[0][0])/dJt;
  
  dJt = sqrt(dJt*dJt);

  G[0][0] = iJt[0][0]*L[0][0] + iJt[0][1]*L[1][0];
  G[0][1] = iJt[0][0]*L[0][1] + iJt[0][1]*L[1][1];
  G[0][2] = iJt[0][0]*L[0][2] + iJt[0][1]*L[1][2];
  G[1][0] = iJt[1][0]*L[0][0] + iJt[1][1]*L[1][0];
  G[1][1] = iJt[1][0]*L[0][1] + iJt[1][1]*L[1][1];
  G[1][2] = iJt[1][0]*L[0][2] + iJt[1][1]*L[1][2];

  GG[0][0] = G[0][0]*G[0][0] + G[1][0]*G[1][0];
  GG[0][1] = G[0][0]*G[0][1] + G[1][0]*G[1][1];
  GG[0][2] = G[0][0]*G[0][2] + G[1][0]*G[1][2];
  GG[1][0] = G[0][1]*G[0][0] + G[1][1]*G[1][0];
  GG[1][1] = G[0][1]*G[0][1] + G[1][1]*G[1][1];
  GG[1][2] = G[0][1]*G[0][2] + G[1][1]*G[1][2];
  GG[2][0] = G[0][2]*G[0][0] + G[1][2]*G[1][0];
  GG[2][1] = G[0][2]*G[0][1] + G[1][2]*G[1][1];
  GG[2][2] = G[0][2]*G[0][2] + G[1][2]*G[1][2];
  
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
    val[0][0] += GG[0][0]*tmp;
    val[0][1] += GG[1][0]*tmp;
    val[0][2] += GG[2][0]*tmp;
    val[1][0] += GG[0][1]*tmp;
    val[1][1] += GG[1][1]*tmp;
    val[1][2] += GG[2][1]*tmp;
    val[2][0] += GG[0][2]*tmp;
    val[2][1] += GG[1][2]*tmp;
    val[2][2] += GG[2][2]*tmp;
  }
  for (ii=0; ii<3; ++ii){
    for (jj=0; jj<3; ++jj){
      val[ii][jj] = val[ii][jj]*dJt;
    }
  }

  return;
}
double phidotphi(double g[3][2], double mua[3], double val[3][3])
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


double boundary_int(double gg[3][2], double bnd_index[3], double ksir[3], \
		    double valr[2][2])
{
  double ksirtmp[2];
  int i, j;
  double g[2][2], val;

  for (i=0; i<2; ++i){
    for (j=0; j<2; ++j){
      valr[i][j]=0;
    }
  }

  if (bnd_index[0]==1 && bnd_index[2]==1){
    ksirtmp[0] = ksir[0];
    ksirtmp[1] = ksir[2];
    g[0][0] = gg[0][0];
    g[0][1] = gg[0][1];
    g[1][0] = gg[2][0];
    g[1][1] = gg[2][1];
    for (i=0; i<2; ++i){
      for (j=0; j<=i; ++j){
	bound2(g, i, j, ksirtmp, &val);
	valr[i][j] = val;
      }
    }
  }
  else if (bnd_index[0]==1 && bnd_index[1]==1){
    ksirtmp[0] = ksir[0];
    ksirtmp[1] = ksir[1];
    g[0][0] = gg[0][0];
    g[0][1] = gg[0][1];
    g[1][0] = gg[1][0];
    g[1][1] = gg[1][1];

    for (i=0; i<2; ++i){
      for (j=0; j<=i; ++j){
	bound2(g, i, j, ksirtmp, &val);
	valr[i][j] = val;
      }
    }
  }
  else if (bnd_index[1]==1 && bnd_index[2]==1){
    ksirtmp[0] = ksir[1];
    ksirtmp[1] = ksir[2];
    g[0][0] = gg[1][0];
    g[0][1] = gg[1][1];
    g[1][0] = gg[2][0];
    g[1][1] = gg[2][1];

    for (i=0; i<2; ++i){
      for (j=0; j<=i; ++j){
	bound2(g, i, j, ksirtmp, &val);
	valr[i][j] = val;
      }
    }
  }
  return;
}
double bound2(double g[2][2], int il, int im, double ksir[2], \
		    double *val)
{
  static double w[2] = {1.0/2.0,1.0/2.0};
  static double ip[2] = {0.21132486540519,0.78867513459481};
  double dJt, S[2], tmp, tmpval;
  int ii;
  double x0 = g[0][0]; double x1 = g[1][0];
  double y0 = g[0][1]; double y1 = g[1][1];
  double k0 = ksir[0]; double k1 = ksir[1];

  dJt = sqrt(((g[1][0]-g[0][0])*(g[1][0]-g[0][0])) +	\
	     ((g[1][1]-g[0][1])*(g[1][1]-g[0][1])));


  tmpval = 0.0;
  for (ii=0; ii<2; ++ii){
    S[0] = 1-ip[ii];
    S[1] = ip[ii]; 
    tmp = k0*S[0] + k1*S[1];
    tmp = tmp*w[ii];
    tmpval += S[il]*S[im]*tmp;
    }
  *val = tmpval*dJt;
}
/* -------- Gate-way to matlab  ------------ */

void mexFunction(int nlhs,
		 mxArray *plhs[],
		 int nrhs,
		 const mxArray *prhs[])
     
{
  double *nodes,*elements,*bndvtx,*mua,*kappa,*ksir,*c,*dt;
  int nodem,noden,elemm,elemn,nzmax;
  double *I,*J,*S1,*S2;
  
  /* Error checking  */
  
  if (nrhs < 7 )
    mexErrMsgTxt(" There is not enough input arguments");
  
  if(nlhs!=4)
    mexErrMsgTxt("This routine requires Three ouput arguments");
  
  nodes=mxGetPr(prhs[0]);          /*  nodes of mesh */
  elements=mxGetPr(prhs[1]);       /*  elements of mesh */
  bndvtx=mxGetPr(prhs[2]);         /*  boundary nodes */
  mua=mxGetPr(prhs[3]);            /*  absorption, nosal */
  kappa=mxGetPr(prhs[4]);          /*  diffusion, nodeal */
  ksir=mxGetPr(prhs[5]);           /*  reflection parameter, nodal  */
  dt=mxGetPr(prhs[6]);            /* dt */
  nodem=mxGetM(prhs[0]);          /*  Number of of nodes */
  noden=mxGetN(prhs[0]);          /*  Number of node freedom */
  elemm=mxGetM(prhs[1]);          /*  Number of elements */
  elemn=mxGetN(prhs[1]);          /*  Number of nodes per elements */
  
  nzmax = elemm*(elemn*(elemn+1));
  /*nzmax = nodem*nodem*0.1;*/
  /*mexPrintf("%d\n",nzmax);*/

  plhs[0]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /* i */
  plhs[1]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /* j */
  plhs[2]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /* s1 */
  plhs[3]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /* s2 */

  I=mxGetPr(plhs[0]);
  J=mxGetPr(plhs[1]);
  S1=mxGetPr(plhs[2]);
  S2=mxGetPr(plhs[3]);
  
  mainloop(nodes,elements,bndvtx,mua,kappa,ksir,nodem,noden,elemm,elemn,dt,I,J,S1,S2);
  
  return;
}
