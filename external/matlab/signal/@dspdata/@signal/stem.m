function varargout = stem(this, varargin)
%STEM   Create a stem plot of the signal.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

h = plotfcn(this, 'stem', varargin{:});

if nargout
    varargout = {h};
end

% [EOF]
