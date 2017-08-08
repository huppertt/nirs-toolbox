function varargout = ellip(this, varargin)
%ELLIP   Elliptic or Cauer digital filter design.
%   H = ELLIP(D) Design an Elliptic digital filter using the specifications
%   in the object D.
%
%   H = ELLIP(D, MATCH) Design a filter and match one band exactly.  MATCH
%   can be either 'passband' 'stopband' or 'both' (default).  This flag is
%   only used when designing minimum order Elliptic filters.

%   Copyright 2005-2009 The MathWorks, Inc.

validate_iir_designmethod(this,'elliptical')

try
    [varargout{1:nargout}] = privdesigngateway(this, 'ellip', ...
        'DesignMode','Interpolator',varargin{:});
catch e
    throw(e);
end
