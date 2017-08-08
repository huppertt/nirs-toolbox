function b = thisisreal(this)
%THISISREAL   True if the object is real.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

Hd = dispatch(this);

b = isreal(Hd);

% [EOF]
