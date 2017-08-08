function this = nyqordntw(varargin)
%NYQORDNTW   Construct a NYQORDNTW object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

this = fspecs.nyqordntw;

this.ResponseType = 'Nyquist with filter order and transition width';

this.setspecs(varargin{:});

% [EOF]
