function varargout = currentfdesigndesignmethods(this,varargin)
%CURRENTFDESIGNDESIGNMETHODS

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

% Return {} if the user requested IIR.
if any(strcmpi(varargin, 'iir'))
    varargout = {{}, false, 'iir'};
    return;
end

% Ask the contained object which FIR design methods are available.
[varargout{1:nargout}] = thisdesignmethods(this, ...
    varargin{:}, 'fir');

% [EOF]
