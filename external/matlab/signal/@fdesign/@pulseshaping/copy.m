function h = copy(this)
%COPY   Copy the designer.

%   Copyright 2008 The MathWorks, Inc.

h = fdesign.pulseshaping;

h.PulseShape = this.PulseShape;
h.PulseShapeObj = copy(this.PulseShapeObj);

% [EOF]
