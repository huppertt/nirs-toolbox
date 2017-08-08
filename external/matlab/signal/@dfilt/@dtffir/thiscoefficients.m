function C = thiscoefficients(Hd)
%THISCOEFFICIENTS Filter coefficients.
%   C = THISCOEFFICIENTS(Hd) returns a cell array of coefficients of
%   discrete-time filter Hd.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.
  
C = {Hd.Numerator};
