function varargout = maxflat(this, varargin)
%MAXFLAT   FIR filter design using the maxflat method
%   MAXFLAT(D) FIR filter design using the maxflat method.

%   Copyright 2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'maxflat', varargin{:});
catch ME
    throwAsCaller(ME);
end

% [EOF]
