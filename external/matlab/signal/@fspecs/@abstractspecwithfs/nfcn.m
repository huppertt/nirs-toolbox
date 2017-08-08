function varargout = nfcn(this, fcn, varargin)
%NFCN   Evaluate a function with normalized frequency set to true.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

normalized = get(this, 'NormalizedFrequency');
normalizefreq(this, true);

[varargout{1:nargout}] = feval(fcn, this, varargin{:});

normalizefreq(this, normalized);

% [EOF]
