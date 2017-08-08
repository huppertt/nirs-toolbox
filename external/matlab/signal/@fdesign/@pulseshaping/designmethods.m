function varargout = designmethods(this, varargin)
%DESIGNMETHODS   Returns a cell of design methods.

%   Copyright 2008 The MathWorks, Inc.

if nargout
    varargout = {designmethods(this.PulseShapeObj, varargin{:})};
else
    designmethods(this.PulseShapeObj, varargin{:})
end
% [EOF]