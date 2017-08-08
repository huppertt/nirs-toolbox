function Hd = dispatch(this)
%DISPATCH   Returns the LWDFILT.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

den = [1, this.AllpassCoefficients];
num = fliplr(den);

Hd = lwdfilt.tf(num,den);

Hd.refNum = num;
Hd.refDen = den;

% [EOF]
