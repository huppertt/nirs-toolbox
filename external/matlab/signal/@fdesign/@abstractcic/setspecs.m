function setspecs(this, M, varargin)
%SETSPECS   Set the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if nargin > 1
    set(this, 'DifferentialDelay', M);
end

abstract_setspecs(this, varargin{:});

% [EOF]
