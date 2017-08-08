function this = nyqmin(varargin)
%NYQMIN   Construct a NYQMIN object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

this = fspecs.nyqmin;

this.ResponseType = 'Minimum-order nyquist';

this.setspecs(varargin{:});

% [EOF]
