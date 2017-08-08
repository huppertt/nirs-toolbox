function this = nyqord(varargin)
%NYQORD   Construct a NYQORD object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

this = fspecs.nyqord;

this.ResponseType = 'Nyquist with filter order';

this.setspecs(varargin{:});

% [EOF]
