function varargout = iirls(this,varargin)
%IIRLS   

%   Author(s): V. Pellissier
%   Copyright 2005-2008 The MathWorks, Inc.


try
    [varargout{1:nargout}] = privdesigngateway(this, 'iirls', varargin{:});
catch e
    throw(e);
end

% [EOF]
