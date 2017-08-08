function Hd = dispatch(this)
%DISPATCH   Return the LWDFILT.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

[b, a] = latc2tf(this.Lattice, this.Ladder);

Hd = lwdfilt.tf(b, a);

[b, a] = latc2tf(this.refLattice, this.refLadder);

Hd.refnum = b;
Hd.refden = a;

% [EOF]
