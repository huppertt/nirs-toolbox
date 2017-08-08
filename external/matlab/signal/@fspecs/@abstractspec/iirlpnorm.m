function varargout = iirlpnorm(this,varargin)
%IIRLPNORM   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

[varargout{1:nargout}] = design(this, 'iirlpnorm', varargin{:});


% [EOF]
