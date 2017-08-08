function ladder = setladder(Hd, ladder)
%SETLADDER Overloaded set on the Ladder property.
  
%   Author: V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

ladder = set_coeffs(Hd, ladder);

% Check data type and store Ladder as reference
set(Hd,'refladder',ladder);

ncoeffs = Hd.ncoeffs;
oldlength=0;
if length(ncoeffs)==2, oldlength = ncoeffs(2); end
newlength = length(ladder);
ncoeffs(2) = newlength;
Hd.ncoeffs = ncoeffs;

% Update the ncoeffs property of the plugin.
set_ncoeffs(Hd.filterquantizer, ncoeffs);

if newlength~=oldlength,
    reset(Hd);
end

% Quantize the coefficients
quantizecoeffs(Hd);

% Hold an empty to not duplicate storage
ladder = [];
