function Hd = todf1tsos(this)
%TODF1TSOS   Convert to a DF1TSOS.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[sos, g] = tf2sos(this.Numerator, this.Denominator);

Hd = dfilt.df1tsos(sos, g);

% [EOF]
