function varargout = iirlpnorm(this,varargin)
%IIRLPNORM   

%   Author(s): V. Pellissier
%   Copyright 2005-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'iirlpnorm', varargin{:});
catch e
    throw(e);
end

% [EOF]
