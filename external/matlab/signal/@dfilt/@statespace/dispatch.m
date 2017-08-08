function Hd = dispatch(this)
%DISPATCH   Return the LWDFILT.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

[b, a] = ss2tf(this.A, this.B, this.C, this.D);

Hd = lwdfilt.tf(b, a);

% [EOF]
