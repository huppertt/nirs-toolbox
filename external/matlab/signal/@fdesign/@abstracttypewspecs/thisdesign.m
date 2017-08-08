function varargout = thisdesign(this, method, varargin)
%THISDESIGN   Design the filter.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

[varargout{1:nargout}] = feval(method, this.CurrentSpecs, varargin{:});

% [EOF]
