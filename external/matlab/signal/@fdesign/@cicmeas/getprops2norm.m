function props2norm = getprops2norm(this)
%GETPROPS2NORM   Get the props2norm.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

props2norm = get(this, { ...
    'Fpass' ...
    'Fstop' ...
    'Fnulls'});

% [EOF]
