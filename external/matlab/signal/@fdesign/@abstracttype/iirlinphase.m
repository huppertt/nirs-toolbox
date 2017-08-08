function varargout = iirlinphase(this,varargin)
%IIRLINPHASE   

%   Author(s): R. Losada
%   Copyright 2005-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'iirlinphase', varargin{:});
catch e
    throw(e);
end



% [EOF]
