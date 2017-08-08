function f = thisisreal(Hd)
%THISISREAL  True for filter with real coefficients.
%   THISISREAL(Hd) returns 1 if filter Hd has real coefficients, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.
  
f = isreal(Hd.sosmatrix) & isreal(Hd.ScaleValues);
