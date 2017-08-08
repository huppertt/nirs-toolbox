function Hd = dispatch(this)
%DISPATCH   Return the lightweight DFILT.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

[b, a] = latc2tf(this.Lattice, 'allpass');

Hd = lwdfilt.tf(b, a);

[b, a] = latc2tf(this.refLattice, 'allpass');

set(Hd, 'refnum', b, 'refden', a);

% [EOF]
