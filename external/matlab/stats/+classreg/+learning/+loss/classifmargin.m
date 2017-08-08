function m = classifmargin(C,Sfit)
%CLASSIFMARGIN Classification margins.
%   M=CLASSIFMARGIN(C,S) returns a column-vector of N classification
%   margins for N rows in C and S. C is an N-by-K logical matrix for N
%   observations and K classes with one true in each row indicating the
%   true class for this observation. S is an N-by-K floating-point matrix
%   of predicted scores. Classification margin is the score for the true
%   class minus the largest score across the false classes.

%   Copyright 2010-2014 The MathWorks, Inc.

[N,K] = size(C);

if K==1
    m = NaN(N,1);
    return;
end

[~,trueC] = max(C,[],2);
m = zeros(N,1);
for k=1:K
    trueK = false(K,1);
    trueK(k) = true;
    idx = trueC==k;
    m(idx) = Sfit(idx,trueK) - max(Sfit(idx,~trueK),[],2);
end
end
