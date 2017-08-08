function g = thisnormalize(Hd)
%THISNORMALIZE   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

g = Hd.refgain;
Hd.refgain = 1;

% [EOF]
