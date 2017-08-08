function Hd = dispatch(this)
%DISPATCH   Return the lwdfilt.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

Hd = lwdfilt.tf(this.Gain);

Hd.refnum = this.refGain;

% [EOF]
