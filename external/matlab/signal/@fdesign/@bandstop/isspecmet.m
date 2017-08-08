function b = isspecmet(this, Hd)
%ISSPECMET   True if the object's specification has been met by the filter.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

b = base_isspecmet(this, Hd, {'Apass1', '<'}, {'Astop', '>'}, {'Apass2', '<'});

% [EOF]
