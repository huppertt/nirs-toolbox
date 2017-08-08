function Hd = dispatch(this)
%DISPATCH   Returns a LWDFILT.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

[b, a] = latc2tf(this.Lattice,'max');

Hd = lwdfilt.tf(b, a);

[b, a] = latc2tf(this.refLattice,'max');

set(Hd, 'refnum', b, 'refden', a);

% [EOF]
