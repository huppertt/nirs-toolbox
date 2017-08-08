function varargout = ifir(this, varargin)
%IFIR    Design an two-stage equiripple filter.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

[varargout{1:nargout}] = design(this, 'ifir', varargin{:});

% [EOF]
