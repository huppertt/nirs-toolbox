function c = thiscoefficients(Hd)
%THISCOEFFICIENTS Filter coefficients.
%   C = THISCOEFFICIENTS(Hd) returns a cell array of coefficients of
%   discrete-time filter Hd.
%
%   See also DFILT.   

%   Author: R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

c = {Hd.sosmatrix, Hd.ScaleValues};
