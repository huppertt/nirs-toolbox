function p = fpval(x,df1,df2)
%FPVAL F distribution p-value function.
%   P = FPVAL(X,V1,V2) returns the upper tail of the F cumulative distribution
%   function with V1 and V2 degrees of freedom at the values in X.  If X is
%   the observed value of an F test statistic, then P is its p-value.
%
%   The size of P is the common size of the input arguments.  A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also FCDF, FINV.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.

%   Copyright 2010 The MathWorks, Inc. 


if nargin < 3, 
    error(message('stats:fpval:TooFewInputs')); 
end

xunder = 1./max(0,x);
xunder(isnan(x)) = NaN;
p = fcdf(xunder,df2,df1);
