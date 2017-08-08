function Hd = firls(this, varargin)
%FIRLS   Design a FIR Least-Squares filter.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

Hd = design(this, 'firls', varargin{:});

% [EOF]
