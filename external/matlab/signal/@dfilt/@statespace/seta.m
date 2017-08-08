function A = seta(Hd, A)
%SETA Overloaded set function on the A property.
  
%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

% Store A so that nstates works properly (needed for states)
set(Hd,'refA',A);
quantizecoeffs(Hd);

% This is not the number of coefficients but this is the number of states
% and the nstates method will use this property.
ncoeffs = Hd.ncoeffs;
oldlength=0;
if ~isempty(ncoeffs), oldlength = ncoeffs(1); end
newlength = size(A,2);
Hd.ncoeffs = newlength;

if newlength~=oldlength,
    reset(Hd);
end


% Hold an empty to not duplicate storage
A = [];

