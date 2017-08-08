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

b = Hd.privNum;
a = Hd.privDen;
nb = length(b);

% TapIndex is a 2 element vector
tapIndexi = Hd.TapIndex;
% In R13, tapIndex was a scalar.
if length(tapIndexi == 1), tapIndexi = [tapIndexi 0]; end

% Extract the numerator and denominator states
ziNum = zi.Numerator;
ziDen = zi.Denominator;

% Returning the numerator and denominator states so that we can get the
% final tap indices and store them back in Hd (since the tapIndex is a
% private property which can only be accessed in a method).
[y,zfNum,zfDen,nBPtrf,dBPtrf] = df1filter(q,b,a,x,ziNum,ziDen,...
    tapIndexi(1),tapIndexi(2));

Hd.TapIndex = [nBPtrf dBPtrf];
zf = filtstates.dfiir(zfNum,zfDen);