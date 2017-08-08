function coeffs = setbeta(Hd,coeffs)
%SETBETA Overloaded set on the Beta property
%
%   This should be a private method.
  
%   Author: V. Pellissier
%   Copyright 1999-2003 The MathWorks, Inc.

set(Hd,'refBeta',coeffs);
quantizecoeffs(Hd);

% Hold an empty to not duplicate storage
coeffs = [];
