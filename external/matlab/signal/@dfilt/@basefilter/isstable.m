function f = isstable(Hb)
%ISSTABLE True if the filter is stable

%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

f = base_is(Hb, 'thisisstable');

% [EOF]
