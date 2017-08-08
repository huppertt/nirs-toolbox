function varargout = butter(this, varargin)
%BUTTER   Butterworth IIR digital filter design.
%   H = BUTTER(D) Design a Butterworth IIR digital filter using the
%   specifications in the object D.
%
%   H = BUTTER(D, MATCH) Design a filter and match one band exactly.  MATCH
%   can be either 'passband' or 'stopband' (default).  This flag is only
%   used when designing minimum order Butterworth filters.

%   Copyright 2005-2009 The MathWorks, Inc.

validate_iir_designmethod(this, 'Butterworth')

try
    [varargout{1:nargout}] = privdesigngateway(this, 'butter',...
        'DesignMode','Decimator',varargin{:});
catch e
    throw(e);
end

% [EOF]
