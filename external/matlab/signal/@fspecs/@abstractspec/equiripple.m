function varargout = equiripple(this, varargin)
%EQUIRIPPLE   Design an equiripple filter.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

[varargout{1:nargout}] = design(this, 'equiripple', varargin{:});

% [EOF]
