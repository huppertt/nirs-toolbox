function proplist = reorderprops(this)
%REORDERPROPS   List of properties to reorder.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

proplist = {'Name','Data',getrangepropname(this)};

% [EOF]
