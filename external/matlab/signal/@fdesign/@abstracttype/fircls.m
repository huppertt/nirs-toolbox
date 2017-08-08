function varargout = fircls(this, varargin)
%FIRCLS   FIR filter design using the constrained least squares method
%   FIRCLS(D) FIR filter design using the constrained least squares method.

%   Copyright 2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'fircls', varargin{:});
catch ME
    throwAsCaller(ME);
end

% [EOF]

