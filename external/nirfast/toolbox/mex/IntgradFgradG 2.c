#include <stdio.h>
#include <math.h>
#include <mex.h>
double Intdd(int i, int j, int k, \
        double g00, double g01, \
        double g10, double g11, \
        double g20, double g21);
/* -------- Heart of the mex file----------- */
void mainloop(double *nodes,
        double *elements,
        double *element_area,
        double *dphir,
        double *dphii,
        double *aphir,
        double *aphii,
        double *Jr,
        double *Ji,
        int nodem,
        int noden,
        int elemm,
        int elemn,
        int isComplex) {
    int ele, i, j, k, indextmp2[3], index[3];
    double tri_vtx[3][2], indextmp, val;
    
    if (isComplex > 0){
        for (ele=0; ele<elemm; ++ele){
            for (i=0; i<elemn; ++i){
                indextmp = *(elements+(ele+(i*elemm)));
                indextmp2[i] = indextmp - 1;
            }
            
            /* sorting index of nodes */
            if (indextmp2[0]<indextmp2[1]) {
                if (indextmp2[1]<indextmp2[2]) {
                    index[0] = indextmp2[0];
                    index[1] = indextmp2[1];
                    index[2] = indextmp2[2];
                }
                else {
                    if (indextmp2[0]<indextmp2[2]) {
                        index[0] = indextmp2[0];
                        index[1] = indextmp2[2];
                        index[2] = indextmp2[1];
                    }
                    else {
                        index[0] = indextmp2[2];
                        index[1] = indextmp2[0];
                        index[2] = indextmp2[1];
                    }
                }
            }
            else {
                if (indextmp2[1]>indextmp2[2]) {
                    index[0] = indextmp2[2];
                    index[1] = indextmp2[1];
                    index[2] = indextmp2[0];
                }
                else {
                    if (indextmp2[0]<indextmp2[2]) {
                        index[0] = indextmp2[1];
                        index[1] = indextmp2[0];
                        index[2] = indextmp2[2];
                    }
                    else {
                        index[0] = indextmp2[1];
                        index[1] = indextmp2[2];
                        index[2] = indextmp2[0];
                    }
                }
            }
            
            
            
            for (i=0; i<elemn; ++i){
                for (j=0; j<noden; ++j){
                    tri_vtx[i][j] = *(nodes+(index[i]+(j*nodem)));
                }
            }
            /*mexPrintf("%d %d %d\n", index[0], index[1], index[2]);*/
            /*mexPrintf("%f %f\n",tri_vtx[0][0], tri_vtx[0][1]);*/
            /*mexPrintf("%f %f\n",tri_vtx[1][0], tri_vtx[1][1]);*/
            /*mexPrintf("%f %f\n",tri_vtx[2][0], tri_vtx[2][1]);*/
            for (i=0; i<3; ++i){
                for (j=0; j<3; ++j){
                    val = Intdd(i, j, j, tri_vtx[0][0], tri_vtx[0][1], \
                            tri_vtx[1][0], tri_vtx[1][1], \
                            tri_vtx[2][0], tri_vtx[2][1]);
                    Jr[index[i]] += ((*(dphir+index[j]) * *(aphir+index[j])) - \
                            (*(dphii+index[j]) * *(aphii+index[j])))*val;
                    Ji[index[i]] += ((*(dphir+index[j]) * *(aphii+index[j])) + \
                            (*(dphii+index[j]) * *(aphir+index[j])))*val;
                    for (k=0; k<=(j-1); ++k){
                        val = Intdd(i, j, k, tri_vtx[0][0], tri_vtx[0][1], \
                                tri_vtx[1][0], tri_vtx[1][1], \
                                tri_vtx[2][0], tri_vtx[2][1]);
                        Jr[index[i]] += (((*(dphir+index[j]) * *(aphir+index[k])) - \
                                (*(dphii+index[j]) * *(aphii+index[k]))) + \
                                ((*(dphir+index[k]) * *(aphir+index[j])) - \
                                (*(dphii+index[k]) * *(aphii+index[j]))))*val;
                        Ji[index[i]] += (((*(dphir+index[j]) * *(aphii+index[k])) + \
                                (*(dphii+index[j]) * *(aphir+index[k]))) + \
                                ((*(dphir+index[k]) * *(aphii+index[j])) + \
                                (*(dphii+index[k]) * *(aphir+index[j]))))*val;
                    }
                }
            }
        }
    } else {
        
        for (ele=0; ele<elemm; ++ele){
            for (i=0; i<elemn; ++i){
                indextmp = *(elements+(ele+(i*elemm)));
                indextmp2[i] = indextmp - 1;
            }
            
            /* sorting index of nodes */
            if (indextmp2[0]<indextmp2[1]) {
                if (indextmp2[1]<indextmp2[2]) {
                    index[0] = indextmp2[0];
                    index[1] = indextmp2[1];
                    index[2] = indextmp2[2];
                }
                else {
                    if (indextmp2[0]<indextmp2[2]) {
                        index[0] = indextmp2[0];
                        index[1] = indextmp2[2];
                        index[2] = indextmp2[1];
                    }
                    else {
                        index[0] = indextmp2[2];
                        index[1] = indextmp2[0];
                        index[2] = indextmp2[1];
                    }
                }
            }
            else {
                if (indextmp2[1]>indextmp2[2]) {
                    index[0] = indextmp2[2];
                    index[1] = indextmp2[1];
                    index[2] = indextmp2[0];
                }
                else {
                    if (indextmp2[0]<indextmp2[2]) {
                        index[0] = indextmp2[1];
                        index[1] = indextmp2[0];
                        index[2] = indextmp2[2];
                    }
                    else {
                        index[0] = indextmp2[1];
                        index[1] = indextmp2[2];
                        index[2] = indextmp2[0];
                    }
                }
            }
            
            
            
            for (i=0; i<elemn; ++i){
                for (j=0; j<noden; ++j){
                    tri_vtx[i][j] = *(nodes+(index[i]+(j*nodem)));
                }
            }
            /*mexPrintf("%d %d %d\n", index[0], index[1], index[2]);*/
            /*mexPrintf("%f %f\n",tri_vtx[0][0], tri_vtx[0][1]);*/
            /*mexPrintf("%f %f\n",tri_vtx[1][0], tri_vtx[1][1]);*/
            /*mexPrintf("%f %f\n",tri_vtx[2][0], tri_vtx[2][1]);*/
            for (i=0; i<3; ++i){
                for (j=0; j<3; ++j){
                    val = Intdd(i, j, j, tri_vtx[0][0], tri_vtx[0][1], \
                            tri_vtx[1][0], tri_vtx[1][1], \
                            tri_vtx[2][0], tri_vtx[2][1]);
                    Jr[index[i]] += (*(dphir+index[j]) * *(aphir+index[j]))*val;
                    for (k=0; k<=(j-1); ++k){
                        val = Intdd(i, j, k, tri_vtx[0][0], tri_vtx[0][1], \
                                tri_vtx[1][0], tri_vtx[1][1], \
                                tri_vtx[2][0], tri_vtx[2][1]);
                        Jr[index[i]] += ((*(dphir+index[j]) * *(aphir+index[k])) + \
                                (*(dphir+index[k]) * *(aphir+index[j])))*val;
                    }
                }
            }
        }
        
    }
    return;
}

double Intdd(int i, int j, int k, double g00, double g01, double g10, \
        double g11, double g20, double g21) {
    double intdd, Jt[2][2], iJt[2][2], dJt, dJti, G[2][3];
    double foo[3][3], L[2][3];
    
    L[0][0] = -1.0;
    L[0][1] = 1.0;
    L[0][2] = 0.0;
    L[1][0] = -1.0;
    L[1][1] = 0.0;
    L[1][2] = 1.0;
    
    Jt[0][0] = L[0][0]*g00 + L[0][1]*g10 + L[0][2]*g20;
    Jt[0][1] = L[0][0]*g01 + L[0][1]*g11 + L[0][2]*g21;
    Jt[1][0] = L[1][0]*g00 + L[1][1]*g10 + L[1][2]*g20;
    Jt[1][1] = L[1][0]*g01 + L[1][1]*g11 + L[1][2]*g21;
    
    dJt = Jt[0][0]*Jt[1][1] - Jt[0][1]*Jt[1][0];
    dJti = 1.0/dJt;
    
    iJt[0][0] = Jt[1][1]*dJti;
    iJt[0][1] = -Jt[0][1]*dJti;
    iJt[1][0] = -Jt[1][0]*dJti;
    iJt[1][1] = Jt[0][0]*dJti;
    
    G[0][0] = iJt[0][0]*L[0][0] + iJt[0][1]*L[1][0];
    G[0][1] = iJt[0][0]*L[0][1] + iJt[0][1]*L[1][1];
    G[0][2] = iJt[0][0]*L[0][2] + iJt[0][1]*L[1][2];
    G[1][0] = iJt[1][0]*L[0][0] + iJt[1][1]*L[1][0];
    G[1][1] = iJt[1][0]*L[0][1] + iJt[1][1]*L[1][1];
    G[1][2] = iJt[1][0]*L[0][2] + iJt[1][1]*L[1][2];
    
    if (dJt < 0)
        dJt = -1.0*dJt;
    dJti = dJt*0.5;
    
    foo[0][0] = (G[0][0]*G[0][0] + G[1][0]*G[1][0])*dJti;
    foo[0][1] = (G[0][0]*G[0][1] + G[1][0]*G[1][1])*dJti;
    foo[0][2] = (G[0][0]*G[0][2] + G[1][0]*G[1][2])*dJti;
    foo[1][0] = (G[0][1]*G[0][0] + G[1][1]*G[1][0])*dJti;
    foo[1][1] = (G[0][1]*G[0][1] + G[1][1]*G[1][1])*dJti;
    foo[1][2] = (G[0][1]*G[0][2] + G[1][1]*G[1][2])*dJti;
    foo[2][0] = (G[0][2]*G[0][0] + G[1][2]*G[1][0])*dJti;
    foo[2][1] = (G[0][2]*G[0][1] + G[1][2]*G[1][1])*dJti;
    foo[2][2] = (G[0][2]*G[0][2] + G[1][2]*G[1][2])*dJti;
    
    if (j>=k)
        intdd = foo[j][k]*0.333333333;
    else
        intdd = foo[k][j]*0.333333333;
    
    return intdd;
}

/* -------- Gate-way to matlab  ------------ */

void mexFunction(int nlhs,
        mxArray *plhs[],
        int nrhs,
        const mxArray *prhs[])
        
{
    double *nodes, *elements, *element_area, *dphir, *dphii, *aphir, *aphii;
    int nodem, noden, elemm, elemn, isComplex;
    double *Jr, *Ji;
    
    /* Error checking  */
    
    if (nrhs < 5 )
        mexErrMsgTxt(" There is not enough input arguments");
    
    if(nlhs!=1)
        mexErrMsgTxt("This routine requires one ouput arguments");
    
    nodes=mxGetPr(prhs[0]);          /*  nodes of mesh */
    elements=mxGetPr(prhs[1]);       /*  elements of mesh */
    element_area=mxGetPr(prhs[2]);   /*  area of elements */
    
    isComplex = mxIsComplex(prhs[3]) + mxIsComplex(prhs[4]);
    
    dphir=mxGetPr(prhs[3]);           /*  direct solution */
    
    if (mxIsComplex(prhs[3])){
        dphii=mxGetPi(prhs[3]);
    }else if (isComplex > 0){							/*  if there is no complex part, make a new array with zeros */ 
        dphii=(double*) mxCalloc( mxGetM(prhs[3]), sizeof(double) );
    }else {
        dphii = 0;
    }
    
    aphir=mxGetPr(prhs[4]);           /*  adjoint solution */
    
    if (mxIsComplex(prhs[4])){
        aphii=mxGetPi(prhs[4]);
    }else if (isComplex > 0){							/*  if there is no complex part, make a new array with zeros */ 
        aphii=(double*) mxCalloc( mxGetM(prhs[4]), sizeof(double) );
    } else {
        aphii = 0;
    }
    
    nodem=mxGetM(prhs[0]);          /*  Number of rows of nodes */
    noden=mxGetN(prhs[0]);          /*  Number of cols of nodes */
    elemm=mxGetM(prhs[1]);           /*  Number of rows of elements */
    elemn=mxGetN(prhs[1]);           /*  Number of cols of elements */
    
    
    if (isComplex > 0){
    plhs[0]=mxCreateDoubleMatrix(nodem, 1, mxCOMPLEX); /* vector to return to Matlab */
    } else {
        plhs[0]=mxCreateDoubleMatrix(nodem, 1, mxREAL);
    }
    
    Jr=mxGetPr(plhs[0]);
    if (isComplex > 0){
    Ji=mxGetPi(plhs[0]);
    }else{
        Ji=0;
    }
    
    
    
    mainloop(nodes, elements, element_area, dphir, dphii, aphir, aphii, Jr, Ji, nodem, noden, elemm, elemn, isComplex);
    
    return;
}
