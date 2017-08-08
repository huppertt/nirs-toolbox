function varargout = cheby1(this, varargin)
%CHEBY1   Chebyshev Type I digital filter design.
%   H = CHEBY1(D) Design a Chebyshev Type I digital filter using the
%   specifications in the object D.
%
%   H = CHEBY1(D, MATCH) Design a filter and match one band exactly.  MATCH
%   can be either 'passband' (default) or 'stopband'.  This flag is only
%   used when designing minimum order Chebyshev filters.

%   Author(s): J. Schickler
%   Copyright 1999-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'cheby1', varargin{:});
catch e
    throw(e);
end


% [EOF]
