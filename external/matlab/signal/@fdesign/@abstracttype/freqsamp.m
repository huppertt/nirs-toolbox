function varargout = freqsamp(this,varargin)
%FREQSAMP   

%   Author(s): V. Pellissier
%   Copyright 2005-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'freqsamp', varargin{:});
catch e
    throw(e);
end


% [EOF]
