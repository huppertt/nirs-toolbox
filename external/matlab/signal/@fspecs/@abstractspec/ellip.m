function Hd = ellip(this, varargin)
%ELLIP   Elliptic or Cauer digital filter design.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

Hd = design(this, 'ellip', varargin{:});

% [EOF]
