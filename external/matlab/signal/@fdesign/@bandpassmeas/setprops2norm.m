function setprops2norm(this, props2norm)
%SETPROPS2NORM   Set the props2norm.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

set(this, {'Fstop1', 'F6dB1', 'F3dB1', 'Fpass1', 'Fpass2', 'F3dB2', ...
    'F6dB2', 'Fstop2'}, props2norm);

% [EOF]
