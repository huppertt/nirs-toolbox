function varargout = autoscale(this,x)
%AUTOSCALE   

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

% Verify that the structure support autoscale
verifyautoscalability(this);

if nargout>0,
    that = copy(this);
else
    that = this;
end

for k=1:length(that.Stage),
    that.Stage(k) = autoscale(that.Stage(k),x);
end

if nargout>0,
    varargout{1} = that;
end

% [EOF]
