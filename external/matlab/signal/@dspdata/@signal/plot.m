function varargout = plot(this, varargin)
%PLOT   Plot the signal.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

h = plotfcn(this, 'line');

if nargout
    varargout = {h};
end

% [EOF]
