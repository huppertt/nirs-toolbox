function props2norm = getprops2norm(this)
%GETPROPS2NORM   Get the props2norm.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% Return the frequency values.
props2norm = get(this, {...
    'Fstop' ...
    'F6dB' ...
    'F3dB' ...
    'Fpass'});

% [EOF]
