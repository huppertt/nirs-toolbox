function Hd = cheby2(this, varargin)
%CHEBY2   Chebyshev Type II digital filter design.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

Hd = design(this, 'cheby2', varargin{:});

% [EOF]
