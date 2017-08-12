#include <mex.h>


void mainloop(int nnodes,int nregions,
    double *region,double *L)
{
    double od[100],ctr[100];
    int i,j,chk;
    
    /*Make the lowest region id = 0*/
    chk = region[0];

    for(i=0;i<nnodes;++i)
    {
        if(region[i] < chk)
            chk = region[i];
    }
    
    if(chk != 0)  /*means the region ids start at 1 possibly?*/
    {
        for(i=0;i<nnodes;++i)
            region[i] = region[i] - chk;
    }
    
    /*Now region ids start at 0*/
        
    for(i=0;i<nregions;++i)
        ctr[i]=0;

    
    /*Get off diagonal values*/
    for(i=0;i<nnodes;++i)
    {
		if (nregions >= region[i])
			ctr[(int)region[i]]++;			
    }
    
    for(i=0;i<nregions;++i)
       od[i] = 1. / ctr[i];
    
    /*Now we have the offdiags - make L*/

    for(i=0;i<nnodes;++i)
    {                                

        for(j=0;j<nnodes;++j)
        {
            if(region[i] == region[j] && i != j)
            {
				if (nregions >= region[i])
					L[nnodes*i+j] = -od[(int)region[i]];
				else
					L[nnodes*i+j] = 0;
            }
            
            if(i == j)
                L[nnodes*i+i] = 1.;
            
            
            
        }
        
    }
 


}


void mexFunction(int nlhs,
mxArray *plhs[],
int nrhs,
const mxArray *prhs[])


{
    int m,n,nregions,nnodes;
    double *L,*region;
    
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    nnodes = mxGetScalar(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(nnodes,nnodes,mxREAL);

    L = mxGetPr(plhs[0]);

    nregions = mxGetScalar(prhs[1]);  /*need to make od vars on the fly*/
    region = mxGetPr(prhs[2]);

    mainloop(nnodes,nregions,region,L);
}
