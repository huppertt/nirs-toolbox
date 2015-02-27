#include <stdio.h>
#include <math.h>
#include <mex.h>
double gradphi(double g[4][3], double kappa[4], double val[4][4]);
double phidotphi(double g[4][3], double mua[4], double val[4][4]);
double boundary_int(double gg[4][3], double bnd_index[4], \
		    double ksir[4], int *ii, int *jj, 
		    double *valr);

/* -------- Heart of the mex file----------- */
void mainloop(double *nodes,
	      double *elements,double *bndvtx,
	      double *mua,double *mus,
	      double *g,double *f1,
	      double *f2,double *f3,
	      double *g1,double *g2,
	      double *g3,double *h1,
	      double *h2,double *h3,
	      double *c,double *omega,
	      int nodem,int noden,
	      int elemm,int elemn,
	      double *I,double *J,
	      double *K1,double *K3,
	      double *K5,double *C,
	      double *C23,double *C49,
	      double *C64,double *C8,
	      double *C16,double *C2_59,
	      double *C2_16,double *C2_49,
	      double *C4_9,double *B,
	      double *F1,double *F2,
	      double *F3,double *G1,
	      double *G2,double *G3,
	      double *H1,double *H2,
	      double *H3,double *X,
	      double *Y)
{
    int ele, i, j, index[4], k,b;
    double indextmp, tri_vtx[4][3];
    double bnd_index[4], kappa1_index[4], kappa3_index[4],kappa5_index[4],kappa7_index[4];
    double mus_index[4], g_index[4];
    double c_index[4], mua_index[4],mua2_59_index[4], mua23_index[4],mua49_index[4];
    double mua2_16_index[4],mua2_49_index[4],mua4_9_index[4];
    double mua64_index[4], mua8_index[4],mua16_index[4];
    double k1_val[4][4],k3_val[4][4], k5_val[4][4];
    double c_val[4][4],c23_val[4][4],c49_val[4][4],c64_val[4][4],c8_val[4][4],c16_val[4][4];
    double f1_index[4],f2_index[4],f3_index[4],g1_index[4],g2_index[4],g3_index[4],c2_59_val[4][4], b_val[4][4];
    double h1_index[4],h2_index[4],h3_index[4],h1_valr, h2_valr,h3_valr;
    double c2_49_val[4][4], c2_16_val[4][4],c4_9_val[4][4];
    double f1_valr, f2_valr,f3_valr, g1_valr, g2_valr,g3_valr;

  k = 0;
  for (ele=0; ele<elemm; ++ele){
    for (i=0; i<elemn; ++i){
      indextmp = *(elements+(ele+(i*elemm)));
      index[i] = indextmp;
    }

    for (i=0; i<elemn; ++i){
       mua_index[i] = *(mua+(index[i]-1));
       mus_index[i] = *(mus+(index[i]-1));
       g_index[i] = *(g+(index[i]-1));
       mua23_index[i] = (2./3.)*mua_index[i];
       mua49_index[i] = (4./9.)*mua_index[i];
       mua64_index[i] = (64./225.)*mua_index[i];
       mua8_index[i]  = (8./15.)*mua_index[i];
       mua16_index[i] = (16./45.)*mua_index[i];
       mua2_59_index[i] = (5./9.)*((mua_index[i]+mus_index[i])-(mus_index[i]*(g_index[i]*g_index[i])));
       mua2_16_index[i] = (16./45.)*((mua_index[i]+mus_index[i])-(mus_index[i]*(g_index[i]*g_index[i])));
       mua2_49_index[i] = (4./9.)*((mua_index[i]+mus_index[i])-(mus_index[i]*(g_index[i]*g_index[i])));   
       mua4_9_index[i] = (9./25.)*((mua_index[i]+mus_index[i])-(mus_index[i]*(g_index[i]*g_index[i]*g_index[i]*g_index[i])));   
       kappa1_index[i] = (1./(3.*((mua_index[i]+mus_index[i])-(mus_index[i]*g_index[i]))));
       kappa3_index[i] = (1./(7.*((mua_index[i]+mus_index[i])-(mus_index[i]*g_index[i]*g_index[i]*g_index[i]))));
       kappa5_index[i] = (1./(11.*((mua_index[i]+mus_index[i])-(mus_index[i]*g_index[i]*g_index[i]*g_index[i]*g_index[i]*g_index[i]))));
       c_index[i] = *omega / *(c+(index[i]-1)) * (-1);
       bnd_index[i] = *(bndvtx+(index[i]-1));
       f1_index[i] = *(f1+(index[i]-1));
       f2_index[i] = *(f2+(index[i]-1));
       f3_index[i] = *(f3+(index[i]-1));
       g1_index[i] = *(g1+(index[i]-1));
       g2_index[i] = *(g2+(index[i]-1));
       g3_index[i] = *(g3+(index[i]-1));
       h1_index[i] = *(h1+(index[i]-1));
       h2_index[i] = *(h2+(index[i]-1));
       h3_index[i] = *(h3+(index[i]-1));

/*mexPrintf("%g %g %g %g\n",mua_index[0],mua_index[1],mua_index[2],mua_index[3]); */


     for (j=0; j<noden; ++j){
	tri_vtx[i][j] = *(nodes+(index[i]-1+(j*nodem)));
     }
    }
       gradphi(tri_vtx,kappa1_index,k1_val);
       gradphi(tri_vtx,kappa3_index,k3_val);
       gradphi(tri_vtx,kappa5_index,k5_val);
       phidotphi(tri_vtx,mua_index,c_val);
       phidotphi(tri_vtx,mua23_index,c23_val);
       phidotphi(tri_vtx,mua49_index,c49_val);
       phidotphi(tri_vtx,mua64_index,c64_val);        
       phidotphi(tri_vtx,mua8_index,c8_val);        
       phidotphi(tri_vtx,mua16_index,c16_val);
       phidotphi(tri_vtx,mua2_59_index,c2_59_val);
       phidotphi(tri_vtx,mua2_16_index,c2_16_val);
       phidotphi(tri_vtx,mua2_49_index,c2_49_val);
       phidotphi(tri_vtx,mua4_9_index,c4_9_val);        
       phidotphi(tri_vtx,c_index,b_val);
       f1_valr = 0.0;
       f2_valr = 0.0;
       f3_valr = 0.0;
       g1_valr = 0.0;
       g2_valr = 0.0;
       g3_valr = 0.0;
       h1_valr = 0.0;
       h2_valr = 0.0;
       h3_valr = 0.0;
       
        


/*    mexPrintf("%g %g %g %g\n %g %g %g %g\n %g %g %g %g\n%g %g %g %g\n", \
	  k_val[0][0], k_val[0][1],k_val[0][2], k_val[0][3], \
	  k_val[1][0], k_val[1][1],k_val[1][2], k_val[1][3], \
	  k_val[2][0], k_val[2][1],k_val[2][2], k_val[2][3], \
	  k_val[3][0], k_val[3][1],k_val[3][2], k_val[3][3]);*/
 /*   mexPrintf("%g %g %g %g\n %g %g %g %g\n %g %g %g %g\n%g %g %g %g\n", \ */
/*     	  c_val[0][0], c_val[0][1],c_val[0][2], c_val[0][3], \ */
/* 	  c_val[1][0], c_val[1][1],c_val[1][2], c_val[1][3], \ */
/* 	  c_val[2][0], c_val[2][1],c_val[2][2], c_val[2][3], \ */
/* 	  c_val[3][0], c_val[3][1],c_val[3][2], c_val[3][3]); */
/*    mexPrintf("%g %g %g %g\n %g %g %g %g\n %g %g %g %g\n%g %g %g %g\n", \
	  b_val[0][0], b_val[0][1],b_val[0][2], b_val[0][3], \
	  b_val[1][0], b_val[1][1],b_val[1][2], b_val[1][3], \
	  b_val[2][0], b_val[2][1],b_val[2][2], b_val[2][3], \
	  b_val[3][0], b_val[3][1],b_val[3][2], b_val[3][3]);*/

    for (i=0; i<4; ++i){
      for (j=0; j<4; ++j){
	if (bnd_index[0]+bnd_index[1]+bnd_index[2]+bnd_index[3] != 0){
	  boundary_int(tri_vtx, bnd_index, f1_index, &i, &j, &f1_valr);
	  boundary_int(tri_vtx, bnd_index, f2_index, &i, &j, &f2_valr);
	  boundary_int(tri_vtx, bnd_index, f3_index, &i, &j, &f3_valr);
	  boundary_int(tri_vtx, bnd_index, g1_index, &i, &j, &g1_valr);
	  boundary_int(tri_vtx, bnd_index, g2_index, &i, &j, &g2_valr);
	  boundary_int(tri_vtx, bnd_index, g3_index, &i, &j, &g3_valr);
	  boundary_int(tri_vtx, bnd_index, h1_index, &i, &j, &h1_valr);
	  boundary_int(tri_vtx, bnd_index, h2_index, &i, &j, &h2_valr);
	  boundary_int(tri_vtx, bnd_index, h3_index, &i, &j, &h3_valr); 
	}
       
      	 I[k] = index[i];
         J[k] = index[j];
         K1[k] = k1_val[j][i];
         K3[k] = k3_val[j][i];
         K5[k] = k5_val[j][i];
         C[k] = c_val[j][i];                
         C23[k] = c23_val[j][i];
         C49[k] = c49_val[j][i];
         C64[k] = c64_val[j][i];
         C8[k] = c8_val[j][i];
         C16[k] = c16_val[j][i];                
         C2_59[k] = c2_59_val[j][i];
         C2_16[k] = c2_16_val[j][i];
         C2_49[k] = c2_49_val[j][i];                
         C4_9[k] = c4_9_val[j][i];
         B[k]= b_val[j][i];
         F1[k] = f1_valr;
         F2[k] = f2_valr;
         F3[k] = f3_valr;
         G1[k] = g1_valr;
         G2[k] = g2_valr;
         G3[k] = g3_valr;
         H1[k] = h1_valr;
         H2[k] = h2_valr;
         H3[k] = h3_valr;

         ++k;
     }
    
    }
  }

/*   k=0; */
/*   for (ele=0; ele<elemm; ++ele){ */
/*     for (i=0; i<elemn; ++i){ */
/*       indextmp = *(elements+(ele+(i*elemm))); */
/*       index[i] = indextmp; */
/*     } */
/*     for (i=0; i<elemn; ++i){ */
/*       bnd_index[i] = *(bndvtx+(index[i]-1)); */
/*       f1_index[i] = *(f1+(index[i]-1)); */
/*       f2_index[i] = *(f2+(index[i]-1)); */
/*       f3_index[i] = *(f3+(index[i]-1)); */
/*       g1_index[i] = *(g1+(index[i]-1)); */
/*       g2_index[i] = *(g2+(index[i]-1)); */
/*       g3_index[i] = *(g3+(index[i]-1)); */
/*       h1_index[i] = *(h1+(index[i]-1)); */
/*       h2_index[i] = *(h2+(index[i]-1)); */
/*       h3_index[i] = *(h3+(index[i]-1)); */

/*       for (j=0; j<noden; ++j){ */
/*               tri_vtx[i][j] = *(nodes+(index[i]-1+(j*nodem))); */
/*       } */
/*     } */
/*     for (i=0; i<4; ++i){ */
/*      for (j=0; j<4; ++j){ */
/*        if (bnd_index[0]+bnd_index[1]+bnd_index[2]+bnd_index[3] != 0){ */
/* 	 boundary_int(tri_vtx, bnd_index, f1_index, &i, &j, &f1_valr); */
/* 	 boundary_int(tri_vtx, bnd_index, f2_index, &i, &j, &f2_valr); */
/* 	 boundary_int(tri_vtx, bnd_index, f3_index, &i, &j, &f3_valr); */
/* 	 boundary_int(tri_vtx, bnd_index, g1_index, &i, &j, &g1_valr); */
/* 	 boundary_int(tri_vtx, bnd_index, g2_index, &i, &j, &g2_valr); */
/* 	 boundary_int(tri_vtx, bnd_index, g3_index, &i, &j, &g3_valr); */
/* 	 boundary_int(tri_vtx, bnd_index, h1_index, &i, &j, &h1_valr); */
/* 	 boundary_int(tri_vtx, bnd_index, h2_index, &i, &j, &h2_valr); */
/* 	 boundary_int(tri_vtx, bnd_index, h3_index, &i, &j, &h3_valr);             */

/* 	  if (bnd_index[0]==1 && bnd_index[1]==1 && bnd_index[2]==1){ */
	    
/*                int index1[3];  */
/*               index1[0] = index[0]; */
/*               index1[1] = index[1]; */
/* 	      index1[2] = index[2]; */
/*               for (i=0; i<3; ++i){ */
/*                   for (j=0; j<=i; ++j){ */
/*                      X[k] = index1[i]; */
/*                      Y[k] = index1[j]; */
/*                      F1[k] = f1_valr; */
/*                      F2[k] = f2_valr; */
/*                      F3[k] = f3_valr; */
/*                      G1[k] = g1_valr; */
/*                      G2[k] = g2_valr; */
/*                      G3[k] = g3_valr; */
/*                      H1[k] = h1_valr; */
/*                      H2[k] = h2_valr; */
/*                      H3[k] = h3_valr;         */
/*                       ++k; */
/*                   } */
/*               } */
/*           } */
/*           else if (bnd_index[0]==1 && bnd_index[1]==1 && bnd_index[3]==1){ */
/*               int index1[3]; */
/*               index1[0] = index[0]; */
/*               index1[1] = index[1]; */
/* 	      index1[2] = index[3]; */
/*               for (i=0; i<3; ++i){ */
/*                   for (j=0; j<=i; ++j){ */
/*                      X[k] = index1[i]; */
/*                      Y[k] = index1[j]; */
/*                      F1[k] = f1_valr; */
/*                      F2[k] = f2_valr; */
/*                      F3[k] = f3_valr; */
/*                      G1[k] = g1_valr; */
/*                      G2[k] = g2_valr; */
/*                      G3[k] = g3_valr; */
/*                      H1[k] = h1_valr; */
/*                      H2[k] = h2_valr; */
/*                      H3[k] = h3_valr;         */
/*                      ++k;  */
/*                   } */
/*               } */
/*           } */
/*           else if (bnd_index[0]==1 && bnd_index[2]==1 && bnd_index[3]==1){ */
/*               int index1[3]; */
/*               index1[0] = index[0]; */
/*               index1[1] = index[2]; */
/* 	      index1[2] = index[3]; */
/*               for (i=0; i<3; ++i){ */
/*                   for (j=0; j<=i; ++j){ */
/*                      X[k] = index1[i]; */
/*                      Y[k] = index1[j]; */
/*                      F1[k] = f1_valr; */
/*                      F2[k] = f2_valr; */
/*                      F3[k] = f3_valr; */
/*                      G1[k] = g1_valr; */
/*                      G2[k] = g2_valr; */
/*                      G3[k] = g3_valr; */
/*                      H1[k] = h1_valr; */
/*                      H2[k] = h2_valr; */
/*                      H3[k] = h3_valr;         */
/*                      ++k;    */
/*                   } */
/*               } */
/*           } */
/* 	  else if (bnd_index[1]==1 && bnd_index[2]==1 && bnd_index[3]==1){ */
/*               int index1[3]; */
/*               index1[0] = index[1]; */
/*               index1[1] = index[2]; */
/* 	      index1[2] = index[3]; */
/*               for (i=0; i<3; ++i){ */
/*                   for (j=0; j<=i; ++j){ */
/*                      X[k] = index1[i]; */
/*                      Y[k] = index1[j]; */
/*                      F1[k] = f1_valr; */
/*                      F2[k] = f2_valr; */
/*                      F3[k] = f3_valr; */
/*                      G1[k] = g1_valr; */
/*                      G2[k] = g2_valr; */
/*                      G3[k] = g3_valr; */
/*                      H1[k] = h1_valr; */
/*                      H2[k] = h2_valr; */
/*                      H3[k] = h3_valr;         */
/*                      ++k;       */
/*                   } */
/*               } */
/*           } */
/*         }  */
     
/*      } */
/*     } */

/*   } */

  return;
}

double gradphi(double g[4][3], double kappa[4], double val[4][4])
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
  double Jt[3][3], dJt, iJt[3][3], S[4], G[3][4], GG[4][4], tmp;
  double x0 = g[0][0]; double x1 = g[1][0];
  double x2 = g[2][0]; double x3 = g[3][0];
  double y0 = g[0][1]; double y1 = g[1][1];
  double y2 = g[2][1]; double y3 = g[3][1];
  double z0 = g[0][2]; double z1 = g[1][2];
  double z2 = g[2][2]; double z3 = g[3][2];
  double k0 = kappa[0]; double k1 = kappa[1];
  double k2 = kappa[2]; double k3 = kappa[3];
  
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
  
  iJt[0][0] = +(Jt[1][1]*Jt[2][2] - Jt[1][2]*Jt[2][1])/dJt;
  iJt[0][1] = -(Jt[0][1]*Jt[2][2] - Jt[0][2]*Jt[2][1])/dJt;
  iJt[0][2] = +(Jt[0][1]*Jt[1][2] - Jt[0][2]*Jt[1][1])/dJt;
  iJt[1][0] = -(Jt[1][0]*Jt[2][2] - Jt[1][2]*Jt[2][0])/dJt;
  iJt[1][1] = +(Jt[0][0]*Jt[2][2] - Jt[0][2]*Jt[2][0])/dJt;
  iJt[1][2] = -(Jt[0][0]*Jt[1][2] - Jt[0][2]*Jt[1][0])/dJt;
  iJt[2][0] = +(Jt[1][0]*Jt[2][1] - Jt[1][1]*Jt[2][0])/dJt;
  iJt[2][1] = -(Jt[0][0]*Jt[2][1] - Jt[0][1]*Jt[2][0])/dJt;
  iJt[2][2] = +(Jt[0][0]*Jt[1][1] - Jt[0][1]*Jt[1][0])/dJt;
  
  dJt = sqrt(dJt*dJt);
  
  G[0][0] = iJt[0][0]*L[0][0] + iJt[0][1]*L[1][0] + iJt[0][2]*L[2][0];
  G[0][1] = iJt[0][0]*L[0][1] + iJt[0][1]*L[1][1] + iJt[0][2]*L[2][1];
  G[0][2] = iJt[0][0]*L[0][2] + iJt[0][1]*L[1][2] + iJt[0][2]*L[2][2];
  G[0][3] = iJt[0][0]*L[0][3] + iJt[0][1]*L[1][3] + iJt[0][2]*L[2][3];
  G[1][0] = iJt[1][0]*L[0][0] + iJt[1][1]*L[1][0] + iJt[1][2]*L[2][0];
  G[1][1] = iJt[1][0]*L[0][1] + iJt[1][1]*L[1][1] + iJt[1][2]*L[2][1];
  G[1][2] = iJt[1][0]*L[0][2] + iJt[1][1]*L[1][2] + iJt[1][2]*L[2][2];
  G[1][3] = iJt[1][0]*L[0][3] + iJt[1][1]*L[1][3] + iJt[1][2]*L[2][3];
  G[2][0] = iJt[2][0]*L[0][0] + iJt[2][1]*L[1][0] + iJt[2][2]*L[2][0];
  G[2][1] = iJt[2][0]*L[0][1] + iJt[2][1]*L[1][1] + iJt[2][2]*L[2][1];
  G[2][2] = iJt[2][0]*L[0][2] + iJt[2][1]*L[1][2] + iJt[2][2]*L[2][2];
  G[2][3] = iJt[2][0]*L[0][3] + iJt[2][1]*L[1][3] + iJt[2][2]*L[2][3];
  
  GG[0][0] = G[0][0]*G[0][0] + G[1][0]*G[1][0] + G[2][0]*G[2][0];
  GG[0][1] = G[0][0]*G[0][1] + G[1][0]*G[1][1] + G[2][0]*G[2][1];
  GG[0][2] = G[0][0]*G[0][2] + G[1][0]*G[1][2] + G[2][0]*G[2][2];
  GG[0][3] = G[0][0]*G[0][3] + G[1][0]*G[1][3] + G[2][0]*G[2][3];
  GG[1][0] = G[0][1]*G[0][0] + G[1][1]*G[1][0] + G[2][1]*G[2][0];
  GG[1][1] = G[0][1]*G[0][1] + G[1][1]*G[1][1] + G[2][1]*G[2][1];
  GG[1][2] = G[0][1]*G[0][2] + G[1][1]*G[1][2] + G[2][1]*G[2][2];
  GG[1][3] = G[0][1]*G[0][3] + G[1][1]*G[1][3] + G[2][1]*G[2][3];
  GG[2][0] = G[0][2]*G[0][0] + G[1][2]*G[1][0] + G[2][2]*G[2][0];
  GG[2][1] = G[0][2]*G[0][1] + G[1][2]*G[1][1] + G[2][2]*G[2][1];
  GG[2][2] = G[0][2]*G[0][2] + G[1][2]*G[1][2] + G[2][2]*G[2][2];
  GG[2][3] = G[0][2]*G[0][3] + G[1][2]*G[1][3] + G[2][2]*G[2][3];
  GG[3][0] = G[0][3]*G[0][0] + G[1][3]*G[1][0] + G[2][3]*G[2][0];
  GG[3][1] = G[0][3]*G[0][1] + G[1][3]*G[1][1] + G[2][3]*G[2][1];
  GG[3][2] = G[0][3]*G[0][2] + G[1][3]*G[1][2] + G[2][3]*G[2][2];
  GG[3][3] = G[0][3]*G[0][3] + G[1][3]*G[1][3] + G[2][3]*G[2][3];
  
  for (ii=0; ii<4; ++ii){
    for (jj=0; jj<4; ++jj){
      val[ii][jj]=0;
    }
  }

  for (ii=0; ii<24; ++ii){
    S[0] = 1-ip[ii][0]-ip[ii][1]-ip[ii][2];
    S[1] = ip[ii][0]; 
    S[2] = ip[ii][1];
    S[3] = ip[ii][2];
    tmp = k0*S[0] + k1*S[1] + k2*S[2] + k3*S[3];
    tmp = tmp*w[ii];
    val[0][0] += GG[0][0]*tmp;
    val[0][1] += GG[1][0]*tmp;
    val[0][2] += GG[2][0]*tmp;
    val[0][3] += GG[3][0]*tmp;
    val[1][0] += GG[0][1]*tmp;
    val[1][1] += GG[1][1]*tmp;
    val[1][2] += GG[2][1]*tmp;
    val[1][3] += GG[3][1]*tmp;
    val[2][0] += GG[0][2]*tmp;
    val[2][1] += GG[1][2]*tmp;
    val[2][2] += GG[2][2]*tmp;
    val[2][3] += GG[3][2]*tmp;
    val[3][0] += GG[0][3]*tmp;
    val[3][1] += GG[1][3]*tmp;
    val[3][2] += GG[2][3]*tmp;
    val[3][3] += GG[3][3]*tmp;
  }
  for (ii=0; ii<4; ++ii){
    for (jj=0; jj<4; ++jj){
      val[ii][jj] = val[ii][jj]*dJt;
    }
  }
  return;
}
double phidotphi(double g[4][3], double mua[4], double val[4][4])
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

double boundary_int(double gg[4][3], double bnd_index[4], double ksir[4], \
		    int *ii, int *jj, \
		    double *valr)
{
  int indextmp1[3];
  double g[3][3], a, b, c, s, side_size, sres, ksitmp_r, ksitmp_i;
  if (bnd_index[0]==1 && bnd_index[1]==1 && bnd_index[2]==1){
    static double sd1_intf1ff[4][4]={{1/10.,1/30.,1/30.,0.0},\
				     {1/30.,1/30.,1/60.,0.0},\
				     {1/30.,1/60.,1/30.,0.0},\
				     {0.0,0.0,0.0,0.0}};
    static double sd1_intf2ff[4][4]={{1/30.,1/30.,1/60.,0.0},\
				     {1/30.,1/10.,1/30.,0.0},\
				     {1/60.,1/30.,1/30.,0.0},\
				     {0.0,0.0,0.0,0.0}};
    static double sd1_intf3ff[4][4]={{1/30.,1/60.,1/30.,0.0},\
				     {1/60.,1/30.,1/30.,0.0},\
				     {1/30.,1/30.,1/10.,0.0},\
				     {0.0,0.0,0.0,0.0}};
    static double sd1_intf4ff[4][4]={{0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0}};
    indextmp1[0] = bnd_index[0];
    indextmp1[1] = bnd_index[1];
    indextmp1[2] = bnd_index[2];

    g[0][0] = gg[0][0];
    g[0][1] = gg[0][1];
    g[0][2] = gg[0][2];
    g[1][0] = gg[1][0];
    g[1][1] = gg[1][1];
    g[1][2] = gg[1][2];
    g[2][0] = gg[2][0];
    g[2][1] = gg[2][1];
    g[2][2] = gg[2][2];

    a = sqrt((g[0][0]-g[1][0])*(g[0][0]-g[1][0]) + \
	     (g[0][1]-g[1][1])*(g[0][1]-g[1][1]) + \
	     (g[0][2]-g[1][2])*(g[0][2]-g[1][2]));
    b = sqrt((g[1][0]-g[2][0])*(g[1][0]-g[2][0]) + \
	     (g[1][1]-g[2][1])*(g[1][1]-g[2][1]) + \
	     (g[1][2]-g[2][2])*(g[1][2]-g[2][2]));
    c = sqrt((g[2][0]-g[0][0])*(g[2][0]-g[0][0]) + \
	     (g[2][1]-g[0][1])*(g[2][1]-g[0][1]) + \
	     (g[2][2]-g[0][2])*(g[2][2]-g[0][2]));
    s = 1/2.*(a+b+c);

    side_size = sqrt(s*(s-a)*(s-b)*(s-c));
    
    sres = (sd1_intf1ff[*ii][*jj] + sd1_intf2ff[*ii][*jj]+ \
	    sd1_intf3ff[*ii][*jj] + sd1_intf4ff[*ii][*jj]);

    ksitmp_r = (ksir[0]+ksir[1]+ksir[2])*1/3.;
    
    *valr = sres*side_size*ksitmp_r;  
  }
  else if (bnd_index[0]==1 && bnd_index[3]==1 && bnd_index[1]==1){ 
    static double sd2_intf1ff[4][4]={{1/10.,1/30.,0.0,1/30.},\
				     {1/30.,1/30.,0.0,1/60.},\
				     {0.0,0.0,0.0,0.0},\
				     {1/30.,1/60.,0.0,1/30.}};
    static double sd2_intf2ff[4][4]={{1/30.,1/30.,0.0,1/60.},\
				     {1/30.,1/10.,0.0,1/30.},\
				     {0.0,0.0,0.0,0.0},\
				     {1/60.,1/30.,0.0,1/30.}};
    static double sd2_intf3ff[4][4]={{0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0}};
    static double sd2_intf4ff[4][4]={{1/30.,1/60.,0.0,1/30.},\
				     {1/60.,1/30.,0.0,1/30.},\
				     {0.0,0.0,0.0,0.0},\
				     {1/30.,1/30.,0.0,1/10.}};
    indextmp1[0] = bnd_index[0];
    indextmp1[1] = bnd_index[3];
    indextmp1[2] = bnd_index[1];

    g[0][0] = gg[0][0];
    g[0][1] = gg[0][1];
    g[0][2] = gg[0][2];
    g[1][0] = gg[3][0];
    g[1][1] = gg[3][1];
    g[1][2] = gg[3][2];
    g[2][0] = gg[1][0];
    g[2][1] = gg[1][1];
    g[2][2] = gg[1][2];

    a = sqrt((g[0][0]-g[1][0])*(g[0][0]-g[1][0]) + \
	     (g[0][1]-g[1][1])*(g[0][1]-g[1][1]) + \
	     (g[0][2]-g[1][2])*(g[0][2]-g[1][2]));
    b = sqrt((g[1][0]-g[2][0])*(g[1][0]-g[2][0]) + \
	     (g[1][1]-g[2][1])*(g[1][1]-g[2][1]) + \
	     (g[1][2]-g[2][2])*(g[1][2]-g[2][2]));
    c = sqrt((g[2][0]-g[0][0])*(g[2][0]-g[0][0]) + \
	     (g[2][1]-g[0][1])*(g[2][1]-g[0][1]) + \
	     (g[2][2]-g[0][2])*(g[2][2]-g[0][2]));
    s = 1/2.*(a+b+c);

    side_size = sqrt(s*(s-a)*(s-b)*(s-c));
    
    sres = (sd2_intf1ff[*ii][*jj] + sd2_intf2ff[*ii][*jj]+ \
	    sd2_intf3ff[*ii][*jj] + sd2_intf4ff[*ii][*jj]);

    ksitmp_r = (ksir[0]+ksir[3]+ksir[1])*1/3.;
    
    *valr = sres*side_size*ksitmp_r;  
  }
  else if (bnd_index[0]==1 && bnd_index[2]==1 && bnd_index[3]==1){ 
    static double sd3_intf1ff[4][4]={{1/10.,0.0,1/30.,1/30.},\
				     {0.0,0.0,0.0,0.0},\
				     {1/30.,0.0,1/30.,1/60.},\
				     {1/30.,0.0,1/60.,1/30.}};
    static double sd3_intf2ff[4][4]={{0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0}};
    static double sd3_intf3ff[4][4]={{1/30.,0.0,1/30.,1/60.},\
				     {0.0,0.0,0.0,0.0},\
				     {1/30.,0.0,1/10.,1/30.},\
				     {1/60.,0.0,1/30.,1/30.}};
    static double sd3_intf4ff[4][4]={{1/30.,0.0,1/60.,1/30.},\
				     {0.0,0.0,0.0,0.0},\
				     {1/60.,0.0,1/30.,1/30.},\
				     {1/30.,0.0,1/30.,1/10.}};
    indextmp1[0] = bnd_index[0];
    indextmp1[1] = bnd_index[2];
    indextmp1[2] = bnd_index[3];

    g[0][0] = gg[0][0];
    g[0][1] = gg[0][1];
    g[0][2] = gg[0][2];
    g[1][0] = gg[2][0];
    g[1][1] = gg[2][1];
    g[1][2] = gg[2][2];
    g[2][0] = gg[3][0];
    g[2][1] = gg[3][1];
    g[2][2] = gg[3][2];

    a = sqrt((g[0][0]-g[1][0])*(g[0][0]-g[1][0]) + \
	     (g[0][1]-g[1][1])*(g[0][1]-g[1][1]) + \
	     (g[0][2]-g[1][2])*(g[0][2]-g[1][2]));
    b = sqrt((g[1][0]-g[2][0])*(g[1][0]-g[2][0]) + \
	     (g[1][1]-g[2][1])*(g[1][1]-g[2][1]) + \
	     (g[1][2]-g[2][2])*(g[1][2]-g[2][2]));
    c = sqrt((g[2][0]-g[0][0])*(g[2][0]-g[0][0]) + \
	     (g[2][1]-g[0][1])*(g[2][1]-g[0][1]) + \
	     (g[2][2]-g[0][2])*(g[2][2]-g[0][2]));
    s = 1/2.*(a+b+c);

    side_size = sqrt(s*(s-a)*(s-b)*(s-c));
    
    sres = (sd3_intf1ff[*ii][*jj] + sd3_intf2ff[*ii][*jj]+ \
	    sd3_intf3ff[*ii][*jj] + sd3_intf4ff[*ii][*jj]);

    ksitmp_r = (ksir[0]+ksir[2]+ksir[3])*1/3.;
    
    *valr = sres*side_size*ksitmp_r;  
  }
  else if (bnd_index[0]==1 && bnd_index[3]==1 && bnd_index[2]==1){ 
    static double sd4_intf1ff[4][4]={{0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0},\
				     {0.0,0.0,0.0,0.0}};
    static double sd4_intf2ff[4][4]={{0.0,0.0,0.0,0.0},\
				     {0.0,1/10.,1/30.,1/30.},\
				     {0.0,1/30.,1/30.,1/60.},\
				     {0.0,1/30.,1/60.,1/30.}};
    static double sd4_intf3ff[4][4]={{0.0,0.0,0.0,0.0},\
				     {0.0,1/30.,1/30.,1/60.},\
				     {0.0,1/30.,1/10.,1/30.},\
				     {0.0,1/60.,1/30.,1/30.}};
    static double sd4_intf4ff[4][4]={{0.0,0.0,0.0,0.0},\
				     {0.0,1/30.,1/60.,1/30.},\
				     {0.0,1/60.,1/30.,1/30.},\
				     {0.0,1/30.,1/30.,1/10.}};
    indextmp1[0] = bnd_index[0];
    indextmp1[1] = bnd_index[3];
    indextmp1[2] = bnd_index[2];

    g[0][0] = gg[0][0];
    g[0][1] = gg[0][1];
    g[0][2] = gg[0][2];
    g[1][0] = gg[3][0];
    g[1][1] = gg[3][1];
    g[1][2] = gg[3][2];
    g[2][0] = gg[2][0];
    g[2][1] = gg[2][1];
    g[2][2] = gg[2][2];
 
    a = sqrt((g[0][0]-g[1][0])*(g[0][0]-g[1][0]) + \
	     (g[0][1]-g[1][1])*(g[0][1]-g[1][1]) + \
	     (g[0][2]-g[1][2])*(g[0][2]-g[1][2]));
    b = sqrt((g[1][0]-g[2][0])*(g[1][0]-g[2][0]) + \
	     (g[1][1]-g[2][1])*(g[1][1]-g[2][1]) + \
	     (g[1][2]-g[2][2])*(g[1][2]-g[2][2]));
    c = sqrt((g[2][0]-g[0][0])*(g[2][0]-g[0][0]) + \
	     (g[2][1]-g[0][1])*(g[2][1]-g[0][1]) + \
	     (g[2][2]-g[0][2])*(g[2][2]-g[0][2]));
    s = 1/2.*(a+b+c);

    side_size = sqrt(s*(s-a)*(s-b)*(s-c));
    
    sres = (sd4_intf1ff[*ii][*jj] + sd4_intf2ff[*ii][*jj]+ \
	    sd4_intf3ff[*ii][*jj] + sd4_intf4ff[*ii][*jj]);

    ksitmp_r = (ksir[0]+ksir[3]+ksir[2])*1/3.;
    
    *valr = sres*side_size*ksitmp_r;  
  }
  return;
}




/* -------- Gate-way to matlab  ------------ */

void mexFunction(int nlhs,
		 mxArray *plhs[],
		 int nrhs,
		 const mxArray *prhs[])
     
{
   double *nodes,*elements,*bndvtx,*mua, *mus, *g;
    double *f1,*f2,*f3,*g1,*g2,*g3,*h1,*h2,*h3,*c,*omega;
    int nodem,noden,elemm,elemn,nzmax;
    double *I,*J,*K1,*K3,*K5, *C23, *C49,*C64,*C8, *C16, *C2_59, *B, *X, *Y;
    double *C, *C2_16, *C2_49,*C4_9;
    double *F1, *F2,*F3, *G1, *G2,*G3, *H1, *H2, *H3; 
  
  /* Error checking  */
  
  if (nrhs < 17 )
    mexErrMsgTxt(" There are not enough input arguments");
  
  if(nlhs!=27)
    mexErrMsgTxt("This routine requires 27 ouput arguments");
  
     nodes=mxGetPr(prhs[0]);          /*  nodes of mesh */
    elements=mxGetPr(prhs[1]);       /*  elements of mesh */
    bndvtx=mxGetPr(prhs[2]);         /*  boundary nodes */
    mua=mxGetPr(prhs[3]);            /*  absorption, nodal */
    mus=mxGetPr(prhs[4]);			   /*  scattering*/
    g=mxGetPr(prhs[5]);	  	  	   /*  anisotropy */
    f1=mxGetPr(prhs[6]);           /*  reflection parameter, nodal  */
    f2=mxGetPr(prhs[7]);           /*  reflection parameter, nodal  */
    f3=mxGetPr(prhs[8]);
    g1=mxGetPr(prhs[9]);           /*  reflection parameter, nodal  */
    g2=mxGetPr(prhs[10]);           /*  reflection parameter, nodal  */
    g3=mxGetPr(prhs[11]);
    h1=mxGetPr(prhs[12]);
    h2=mxGetPr(prhs[13]);
    h3=mxGetPr(prhs[14]);
    c=mxGetPr(prhs[15]);              /*  speed of light in tissue, nodal */
    omega=mxGetPr(prhs[16]);          /*  frequency of excitation */
    
    
    nodem=mxGetM(prhs[0]);          /*  Number of of nodes */
    noden=mxGetN(prhs[0]);          /*  Number of node freedom */
    elemm=mxGetM(prhs[1]);          /*  Number of elements */
    elemn=mxGetN(prhs[1]);          /*  Number of nodes per elements */
    
    nzmax = elemm*(elemn*(elemn+1));
  /*mexPrintf("%d\n",nzmax);*/
    
    plhs[0]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /* i */
    plhs[1]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /* j */
    plhs[2]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /* K1 */
    plhs[3]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /*K3*/
    plhs[4]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /*(2/3)C*/
    plhs[5]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /*(4/9)C*/
    plhs[6]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /* (64/225)C*/
    plhs[7]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /* (8/15)C */
    plhs[8]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /*(16/45)C */
    plhs[9]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /*(5/9)C2 */
    plhs[10]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /*(16/45)C2 */
    plhs[11]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    /*(4/9)C2   */ 
    plhs[12]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[13]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    
    plhs[14]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[15]=mxCreateDoubleMatrix(nzmax,1,mxCOMPLEX);
    plhs[16]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[17]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[18]=mxCreateDoubleMatrix(nzmax,1,mxREAL);    
    plhs[19]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[20]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[21]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[22]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[23]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[24]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[25]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[26]=mxCreateDoubleMatrix(nzmax,1,mxREAL);


    I=mxGetPr(plhs[0]);
    J=mxGetPr(plhs[1]);
    K1=mxGetPr(plhs[2]);
    K3=mxGetPr(plhs[3]);
    K5=mxGetPr(plhs[4]);
    C=mxGetPr(plhs[5]);
    C23=mxGetPr(plhs[6]);
    C49=mxGetPr(plhs[7]);
    C64=mxGetPr(plhs[8]);
    C8=mxGetPr(plhs[9]);
    C16=mxGetPr(plhs[10]);
    C2_59= mxGetPr(plhs[11]);
    C2_16=mxGetPr(plhs[12]);
    C2_49=mxGetPr(plhs[13]);
    C4_9=mxGetPr(plhs[14]);
    B=mxGetPi(plhs[15]);
    F1=mxGetPr(plhs[16]);
    F2=mxGetPr(plhs[17]);
    F3=mxGetPr(plhs[18]);    
    G1=mxGetPr(plhs[19]);
    G2=mxGetPr(plhs[20]);
    G3=mxGetPr(plhs[21]);
    H1=mxGetPr(plhs[22]);
    H2=mxGetPr(plhs[23]);
    H3=mxGetPr(plhs[24]);    
    X=mxGetPr(plhs[25]);
    Y=mxGetPr(plhs[26]);
    
  
    mainloop(nodes,elements,bndvtx,mua,mus,g,f1,f2,f3,g1,g2,g3,h1,h2,h3,c,omega,nodem, \
    noden,elemm,elemn,I,J,K1,K3,K5,C,C23,C49,C64,C8,C16,C2_59,C2_16,C2_49,C4_9, \
    B,F1,F2,F3,G1,G2,G3,H1,H2,H3,X,Y);
    
  return;
}
