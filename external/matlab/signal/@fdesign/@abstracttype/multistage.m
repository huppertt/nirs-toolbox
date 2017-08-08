function varargout = multistage(this, varargin)
%MULTISTAGE   Design a multistage FIR filter using the equiripple method.
%   MULTISTAGE(D) designs a multistage equiripple filter using the specifications
%   in the object D.

%   Author(s): R. Losada
%   Copyright 2005-2008 The MathWorks, Inc.

try
    [varargout{1:nargout}] = privdesigngateway(this, 'multistage', varargin{:});
catch e
    throw(e);
end

% [EOF]
