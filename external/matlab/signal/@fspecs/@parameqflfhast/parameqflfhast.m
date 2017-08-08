function this = parameqastflfh(varargin)
%PARAMEQFLFHAST   Construct a PARAMEQFLFHAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.parameqflfhast;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
