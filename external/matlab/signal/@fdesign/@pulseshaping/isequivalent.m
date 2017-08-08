function b = isequivalent(this, htest)
%ISEQUIVALENT   True if the object is equivalent.

%   Copyright 2008 The MathWorks, Inc.

if isa(htest, 'fdesign.pulseshaping'),
    b = isequivalent(this.PulseShapeObj, htest.PulseShapeObj);
else
    b = false;
end

% [EOF]
