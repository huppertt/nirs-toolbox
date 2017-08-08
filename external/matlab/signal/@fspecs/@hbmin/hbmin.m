function this = hbmin(varargin)
%HBMIN   Construct a HBMIN object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

this = fspecs.hbmin;

this.ResponseType = 'Minimum-order halfband';

this.setspecs(varargin{:});

% [EOF]
