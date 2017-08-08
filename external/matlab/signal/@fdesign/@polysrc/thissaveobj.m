function s = thissaveobj(this)
%THISSAVEOBJ Save this object.

%   Copyright 2007 The MathWorks, Inc.

s.InterpolationFactor = get(this, 'InterpolationFactor');
s.DecimationFactor = get(this, 'DecimationFactor');

% [EOF]
