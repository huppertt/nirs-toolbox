function varargout = lagrange(this,varargin)
%LAGRANGE   

%   Author(s): V. Pellissier
%   Copyright 2005-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'lagrange', varargin{:});
catch e
    throw(e);
end

% [EOF]
