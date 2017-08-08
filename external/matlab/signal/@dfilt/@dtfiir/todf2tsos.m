function Hd = todf2tsos(this)
%TODF2TSOS   Convert to a DF2TSOS.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[sos, g] = tf2sos(this.Numerator, this.Denominator);

Hd = dfilt.df2tsos(sos, g);

% [EOF]
