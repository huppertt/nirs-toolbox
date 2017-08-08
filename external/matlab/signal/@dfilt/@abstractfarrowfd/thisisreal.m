function b = thisisreal(this)
%THISISREAL   True if the object is real.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

b = isreal(this.Coefficients);

% [EOF]
