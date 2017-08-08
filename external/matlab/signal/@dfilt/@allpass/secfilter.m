function [y,zf] = secfilter(Hd,x,zi)
%SECFILTER Filter this section.
%   [Y,Zf] = SECFILTER(Hd,X,ZI) filters this section.  This function is only
%   intended to be called from DFILT/FILTER.  The initial conditions have
%   already been padded for the C++ implementation.
%
%   See also DFILT.    

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

  
q = Hd.filterquantizer;
a = Hd.privallpasscoeffs;
[y,zf] = allpassfilter(q,a,x,zi);



% [EOF]
