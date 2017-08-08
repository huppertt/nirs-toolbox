function varargout = validate(this)
%VALIDATE Validate the specifications.

%   Copyright 2010 The MathWorks, Inc.

% Pass validation to the response type.
[varargout{1:nargout}] = validate(this.responseTypeSpecs, this);

% [EOF]
