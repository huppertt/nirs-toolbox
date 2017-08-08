function help(this, designmethod)
%HELP   Provide help for the specified design method.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if nargin < 2
    help('fdesign');
elseif isdesignmethod(this, designmethod)
    help(this.CurrentSpecs, designmethod);
else
    error(message('signal:fdesign:abstracttypewspecs:help:invalidDesignMethod', designmethod));
end

% [EOF]
