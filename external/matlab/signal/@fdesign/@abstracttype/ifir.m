function varargout = ifir(this, varargin)
%IFIR   Design a two-stage FIR filter using the IFIR method.
%   IFIR(D) designs a two-stage equiripple filter using the specifications
%   in the object D.

%   Author(s): R. Losada
%   Copyright 2005-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'ifir', varargin{:});
catch e
    throw(e);
end

% [EOF]
