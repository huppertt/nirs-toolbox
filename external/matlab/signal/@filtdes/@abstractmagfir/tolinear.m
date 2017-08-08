function newval = tolinear(h,val,passOrStop)
%TOLINEAR Convert dB value to linear.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

switch passOrStop,
    
case 'pass',
    newval = (10^(val/20) - 1)/(10^(val/20) + 1);
    
case 'stop',
    newval = 10^(-val/20);
end
