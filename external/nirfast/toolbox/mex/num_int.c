#include <stdio.h>
#include <math.h>
#include <mex.h>


/* -------- Heart of the mex file----------- */
void mainloop(int nnodes,
double *nm,
double *n0,
double *Ca,
double *Rn,
double *I,
double *R,
double *bndvtx)
{
    int N, NN, j, w, k,l,b;
    double m, Ca_tmp;
    
    /*      initialise k*/
    k = 0;
    
    for (N=0; N<15; ++N){
        /* initialise I*/
        for (NN=0; NN<nnodes; ++NN){
            *(I+NN) = 0;
        }
        
        
        for (j=0; j<nnodes; ++j){
            b= *(bndvtx+j);
         /*   mexPrintf("%g\n",b);*/
            if (b==1){
                for(w=0; w<1999; ++w){
                    m= (w+1)/2000.0;
                    
                    Ca_tmp = *(Ca+j);
                    if (m < Ca_tmp){
                        R[j]=1.0;
                    }
                    else {
                        
                        R[j]=0.5*pow((((nm[j]*sqrt(1-pow(nm[j],2)*  \
                        (1-pow(m,2))))-(n0[j]*m))/((nm[j]*sqrt(1-pow(nm[j],2)*\
                        (1-pow(m,2))))+(n0[j]*m))),2)+ 0.5*pow((((nm[j]*m)-  \
                        (n0[j]*sqrt(1-pow(nm[j],2)*(1-pow(m,2)))))/((nm[j]*m)+ \
                        (n0[j]*sqrt(1-pow(nm[j],2)*(1-pow(m,2)))))),2);
                        
                    }
                    
                    I[j]=I[j]+R[j]*pow(m,N)*0.0005;
                }
                
            }
            else{
                I[j]=0;
            }
        }
        
        for (l=0; l<nnodes; ++l){
            
            Rn[k]=I[l];
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
    double *Rn, *nm, *n0,*Ca, *I, *R, *bndvtx;
    int nnodes, nzmax;
    
    
  /* Error checking  */
    
    if (nrhs < 4 )
        mexErrMsgTxt(" There are not enough input arguments");
    
    if(nlhs!=1)
        mexErrMsgTxt("This routine requires one ouput argument");
    
/*  nodes=mxGetPr(prhs[0]);         */
    nm=mxGetPr(prhs[0]);
    n0=mxGetPr(prhs[1]);
    Ca=mxGetPr(prhs[2]);
    bndvtx=mxGetPr(prhs[3]);
    
    nnodes=mxGetM(prhs[0]);          /*  Number of of nodes */
    
    /*   I = mxCreateDoubleMatrix(nnodes,1,mxREAL);*/
    /*   R = mxCreateDoubleMatrix(nnodes,1,mxREAL);*/
    
    nzmax=nnodes*15;
    
    plhs[0]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(nnodes,1,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(nnodes,1,mxREAL);
    Rn=mxGetPr(plhs[0]);
    I=mxGetPr(plhs[1]);
    R=mxGetPr(plhs[2]);
    
    
    mainloop(nnodes,nm,n0,Ca,Rn,I,R, bndvtx);
    
    return;
}

