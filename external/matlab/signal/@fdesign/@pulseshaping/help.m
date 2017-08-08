function help(this, designmethod)
%HELP   Provide help for the specified design method.

%   Copyright 2008 The MathWorks, Inc.

if nargin < 2
    help('fdesign');
elseif isdesignmethod(this.PulseShapeObj, designmethod)
    help(this.PulseShapeObj, designmethod);
else
    error(message('signal:fdesign:pulseshaping:help:invalidDesignMethod', designmethod));
end

% [EOF]
