function varargout = iirlinphase(this,varargin)
%IIRLINPHASE   

%   Copyright 2005-2009 The MathWorks, Inc.

validate_iir_designmethod(this,'IIR linear phase')

try
    [varargout{1:nargout}] = privdesigngateway(this, 'iirlinphase', ...
        'DesignMode','Interpolator',varargin{:});
catch e
    throw(e);
end



% [EOF]
