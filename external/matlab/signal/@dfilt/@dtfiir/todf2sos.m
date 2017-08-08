function Hd = todf2sos(this)
%TODF2SOS   Convert to a DF2SOS.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[sos, g] = tf2sos(this.Numerator, this.Denominator);

Hd = dfilt.df2sos(sos, g);

% [EOF]
