function varargout = freqsamp(this,varargin)
%FREQSAMP   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

[varargout{1:nargout}] = design(this, 'freqsamp', varargin{:});


% [EOF]
