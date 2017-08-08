function Hd = dispatch(this)
%DISPATCH   Dispatch to the lwdfilt object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

Hd = lwdfilt.asymfir(this.Numerator);
Hd.refnum = this.refnum;

% [EOF]
