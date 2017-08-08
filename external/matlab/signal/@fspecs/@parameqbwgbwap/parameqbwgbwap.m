function this = parameqbwgbwap(varargin)
%PARAMEQBWGBWAP   Construct a PARAMEQBWGBWAP object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.parameqbwgbwap;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
