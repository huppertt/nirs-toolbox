function str = gettitlestr(d)
%GETTITLESTR Returns the title for GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% This should be private

str = sprintf('%% %s filter designed using the %s function.', ...
    get(d, 'Tag'), ...
    upper(designfunction(d)));

% [EOF]
