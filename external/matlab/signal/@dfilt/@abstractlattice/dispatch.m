function Hd = dispatch(this)
%DISPATCH   Return the LWDFILT.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

Hd = lwdfilt.tf(latc2tf(this.Lattice));

Hd.refNum = latc2tf(this.refLattice, 'fir');

% [EOF]
