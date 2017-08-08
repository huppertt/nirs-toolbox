function varargout = firmultirate(this, method, varargin)
%FIRMULTIRATE   Perform the design of a multirate.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% This is a separate function so that it can be overloaded in case we need
% to force an L*PL order for minimum order cases.
[varargout{1:nargout}] = feval(method, this, varargin{:});

% [EOF]
