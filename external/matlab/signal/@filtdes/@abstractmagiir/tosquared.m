function newval = tolinear(h,val,notused)
%TOLINEAR Convert dB value to linear.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


newval = 1/(10^(val/10));
