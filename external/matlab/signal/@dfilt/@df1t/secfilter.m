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

% Extract the numerator and denominator states
ziNum = zi.Numerator;
ziDen = zi.Denominator;

[y,zfNum,zfDen] = df1tfilter(q,b,a,x,ziNum,ziDen);
zf = filtstates.dfiir(zfNum,zfDen);
