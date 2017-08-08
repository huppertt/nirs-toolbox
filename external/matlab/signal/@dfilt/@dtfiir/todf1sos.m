function Hd = todf1sos(this)
%TODF1SOS   Convert to a DF1SOS.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[sos, g] = tf2sos(this.Numerator, this.Denominator);

Hd = dfilt.df1sos(sos, g);

% [EOF]
