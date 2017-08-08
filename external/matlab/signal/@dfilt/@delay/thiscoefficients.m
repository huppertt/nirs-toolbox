function c = thiscoefficients(this)
%THISCOEFFICIENTS Filter coefficients.
%   C = THISCOEFFICIENTS(Hd) returns a cell array of coefficients of
%   discrete-time filter Hd.
%
  
%   Copyright 1988-2005 The MathWorks, Inc.

Hd = dispatch(this);
c = {Hd.Numerator};
