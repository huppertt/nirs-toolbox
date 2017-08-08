function Hd = dispatch(this)
%DISPATCH   Dispatch to a light weight dfilt.

%   Author(s): J. Schickler
%   Copyright 2004-2008 The MathWorks, Inc.

checksv(this);

Hd = lwdfilt.sos(this.SOSMatrix, this.ScaleValues);

Hd.refSOSMatrix   = this.refsosMatrix;
Hd.refScaleValues = this.refScaleValues;

% [EOF]
