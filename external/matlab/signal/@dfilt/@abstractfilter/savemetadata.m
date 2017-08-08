function s = savemetadata(this)
%SAVEMETADATA   Save the metadata.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

s.fdesign      = getfdesign(this);
s.fmethod      = getfmethod(this);
s.measurements = get(this, 'privMeasurements');

% [EOF]
