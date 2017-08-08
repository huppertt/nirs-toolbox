function varargout = butter(this, varargin)
%BUTTER   Butterworth IIR digital filter design.
%   H = BUTTER(D) Design a Butterworth IIR digital filter using the
%   specifications in the object D.
%
%   H = BUTTER(D, MATCH) Design a filter and match one band exactly.  MATCH
%   can be either 'passband' or 'stopband' (default).  This flag is only
%   used when designing minimum order Butterworth filters.

%   Author(s): J. Schickler
%   Copyright 1999-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'butter', varargin{:});
catch e
    throw(e);
end

% [EOF]
