function coeffs = setallpass1(Hd,coeffs)
%SETALLPASS1 Overloaded set on the Allpass1 property
%
%   This should be a private method.
  
%   Author: V. Pellissier
%   Copyright 1999-2003 The MathWorks, Inc.

ncoeffs = Hd.ncoeffs;
oldlength=0;
if ~isempty(ncoeffs), oldlength = ncoeffs(1); end
newlength = length(coeffs);
ncoeffs(1) = newlength;
Hd.ncoeffs = ncoeffs;

if newlength~=oldlength,
    reset(Hd);
end

set(Hd,'refAllpass1',coeffs);
quantizecoeffs(Hd);

% Hold an empty to not duplicate storage
coeffs = [];
