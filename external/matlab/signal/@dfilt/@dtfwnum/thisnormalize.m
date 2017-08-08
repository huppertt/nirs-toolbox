function g = thisnormalize(Hd)
%THISNORMALIZE   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

num = Hd.refnum;
g = max(abs(num));
Hd.refnum= num/g;

% [EOF]
