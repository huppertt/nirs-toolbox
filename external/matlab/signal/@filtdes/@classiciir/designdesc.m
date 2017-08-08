function str = designdesc(d)
%DESIGNDESC   Returns the design comment.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

str = sprintf('%% Construct an FDESIGN object and call its %s method.', ...
    upper(designfunction(d)));

% [EOF]
