function varargout = multistage(this, varargin)
%MULTISTAGE    Design a multistage equiripple filter.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

[varargout{1:nargout}] = design(this, 'multistage', varargin{:});

% [EOF]
