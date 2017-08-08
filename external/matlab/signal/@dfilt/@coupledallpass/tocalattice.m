function Hd = tocalattice(this)
%TOCALATTICE   Convert to the non-pc calattice.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

Hd = dfilt.calattice(this.Allpass1, this.Allpass2, this.Beta);

% [EOF]
