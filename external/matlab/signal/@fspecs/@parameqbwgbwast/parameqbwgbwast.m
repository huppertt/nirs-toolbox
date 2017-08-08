function this = parameqbwgbwast(varargin)
%PARAMEQBWGBWAST   Construct a PARAMEQBWGBWAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.parameqbwgbwast;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
