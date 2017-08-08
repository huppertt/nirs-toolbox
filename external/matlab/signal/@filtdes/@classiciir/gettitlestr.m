function str = gettitlestr(this)
%GETTITLESTR   PreGet function for the 'titlestr' property.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

str = sprintf('%% %s %s filter designed using FDESIGN.%s.', ...
    get(this, 'Tag'), ...
    get(this, 'ResponseType'), ...
    upper(get(this, 'ResponseType')));

% [EOF]
