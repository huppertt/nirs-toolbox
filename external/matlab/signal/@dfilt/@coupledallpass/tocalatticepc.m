function Hd = tocalatticepc(this)
%TOCALATTICEPC   Convert to the pc calattice.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

Hd = dfilt.calatticepc(this.Allpass1, this.Allpass2, this.Beta);

% [EOF]
