function [y,zf] = secfilter(Hd,x,zi)
%SECFILTER Filter this section.
%   [Y,Zf] = SECFILTER(Hd,X,ZI) filters this section.  This function is only
%   intended to be called from DFILT/FILTER.  The initial conditions have
%   already been padded for the C++ implementation.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.
  
q = Hd.filterquantizer;

% TapIndex is a 2 element vector; Since the DF2 only has one circular
% buffer, we will be updating the first element of the vector.
tapIndexi = Hd.TapIndex;

% Build the state vector
nz = max(Hd.ncoeffs(1),Hd.ncoeffs(2));
zi = zi(1:nz,:);

b = Hd.privNum;
a = Hd.privDen;

% Returning the numerator and denominator states so that we can get the
% final tap index and store it back in Hd (since the tapIndex is a
% private property which can only be accessed in a method).
[y,zf,tapidxf] = df2filter(q,b,a,x,zi,tapIndexi(1));

if length(tapIndexi) < 2
    tapIndexi(2) = 0;
end

Hd.TapIndex = [tapidxf tapIndexi(2)];
