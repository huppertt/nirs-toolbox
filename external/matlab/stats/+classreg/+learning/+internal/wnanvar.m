function V = wnanvar(X,W,bias)
%WNANVAR Weighted variance ignoring NaN's.
%   V=WNANVAR(X,W,0) returns max likelihood estimate of variance along 1st
%   dimension of X weighted by vector W.
%
%   V=WNANVAR(X,W,1) returns unbiased estimate of variance along 1st
%   dimension of X weighted by vector W.

%   Copyright 2011 The MathWorks, Inc.


W(isnan(W)) = 0;
W = W(:);

M = classreg.learning.internal.wnanmean(X,W);

X = bsxfun(@times,bsxfun(@minus,X,M).^2,W);
tfnan = isnan(X);
X(tfnan) = 0;

Wcol = sum(bsxfun(@times,~tfnan,W),1);

V = sum(X,1) ./ Wcol;

if size(X,1)>1 && bias
    Wcol2 = sum(bsxfun(@times,~tfnan,W.^2),1);
    V = sum(X,1) ./ (Wcol-Wcol2./Wcol);
end
end
