function Hd = butter(this, varargin)
%BUTTER   Butterworth digital filter design.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

Hd = design(this, 'butter', varargin{:});

% [EOF]
