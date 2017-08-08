function varargout = fircls(this, varargin)
%FIRCLS   Design a constrained least-squares filter.

%   Copyright 2008 The MathWorks, Inc.

[varargout{1:nargout}] = design(this, 'fircls', varargin{:});