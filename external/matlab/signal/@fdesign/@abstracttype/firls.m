function varargout = firls(this, varargin)
%FIRLS   Design a least-squares filter.   

%   Author(s): J. Schickler
%   Copyright 1999-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'firls', varargin{:});
catch e
    throw(e);
end

% [EOF]
