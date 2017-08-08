function c = coefficientvariables(Hb)
%COEFFICIENTVARIABLES Coefficient variables.

%   This should be a private method.

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.

Hd = dispatch(Hb);
c = coefficientvariables(Hd);

% [EOF]
