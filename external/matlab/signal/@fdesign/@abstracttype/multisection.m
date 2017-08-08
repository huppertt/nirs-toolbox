function varargout = multisection(this, varargin)
%MULTISECTION   

%   Author(s): J. Schickler
%   Copyright 2005-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'multisection', varargin{:});
catch e
    throw(e);
end

% [EOF]
