function varargout = window(this, varargin)
%WINDOW   FIR filter design using the window method.
%   WINDOW(D) FIR filter design using the window method.

%   Author(s): J. Schickler
%   Copyright 1999-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'window', varargin{:});
catch e
    throw(e);
end

% [EOF]
