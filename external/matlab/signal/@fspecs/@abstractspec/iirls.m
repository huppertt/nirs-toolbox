function varargout = iirls(this,varargin)
%IIRLS   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

[varargout{1:nargout}] = design(this, 'iirls', varargin{:});


% [EOF]
