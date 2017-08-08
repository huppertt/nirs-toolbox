function str = gettitlestr(d)
%GETTITLESTR Returns the title for GENMCODE

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% This should be private

str = sprintf('%% IIR Notching filter designed using the %s function.', ...
    upper(designfunction(d)));
    
% [EOF]
