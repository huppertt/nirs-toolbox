function str = designdesc(d)
%DESIGNDESC Returns the comment that precedes the design call.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

str = sprintf('%% Calculate the coefficients using the %s function.', ...
    upper(designfunction(d)));

% [EOF]
