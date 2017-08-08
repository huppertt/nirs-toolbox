function sglstates = single(h)
%SINGLE   Convert a FILTSTATES.DFIIR object to a single vector.

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.

sglstates = [single(h.Numerator); single(h.Denominator)];

% [EOF]
