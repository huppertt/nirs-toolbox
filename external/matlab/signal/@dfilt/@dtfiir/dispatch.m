function Hd = dispatch(this)
%DISPATCH   Return the LWDFILT.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

Hd = lwdfilt.tf(this.Numerator, this.Denominator);

set(Hd, 'refnum', this.refnum, ...
    'refden', this.refden);

% [EOF]
