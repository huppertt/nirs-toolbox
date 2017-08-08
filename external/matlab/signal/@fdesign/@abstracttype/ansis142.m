function varargout = ansis142(this,varargin)
%ANSIS142   

%   Copyright 2009 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'ansis142', varargin{:});
catch e
    throw(e);
end

% [EOF]
