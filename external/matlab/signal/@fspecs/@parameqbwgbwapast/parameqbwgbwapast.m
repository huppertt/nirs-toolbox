function this = parameqbwgbwapast(varargin)
%PARAMEQBWGBWAPAST   Construct a PARAMEQBWGBWAPAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.parameqbwgbwapast;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
