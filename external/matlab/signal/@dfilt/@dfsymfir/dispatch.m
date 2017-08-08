function Hd = dispatch(this)
%DISPATCH   Dispatch to the lwdfilt object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

Hd = lwdfilt.symfir(this.Numerator);
Hd.refnum = this.refnum;

% [EOF]
