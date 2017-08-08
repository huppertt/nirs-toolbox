function Hd = dispatch(this)
%DISPATCH   Returns the LWDFILT.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

Hd = lwdfilt.tf(this.Numerator);

Hd.refNum = this.refNum;

% [EOF]
