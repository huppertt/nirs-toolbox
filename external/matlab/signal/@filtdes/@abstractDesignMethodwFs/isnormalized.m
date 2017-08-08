function boolflag = isnormalized(d)
%ISNORMALIZED Returns true if the filtdes is normalized.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

freqopt = set(d, 'freqUnits');
boolflag = strcmp(get(d, 'freqUnits'), freqopt{1});

% [EOF]
