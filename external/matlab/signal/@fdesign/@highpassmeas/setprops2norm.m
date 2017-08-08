function setprops2norm(this, props2norm)
%SETPROPS2NORM   Set the props2norm.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% Set the frequency values.
set(this, {'Fstop', 'F6dB', 'F3dB', 'Fpass'}, props2norm);

% [EOF]
