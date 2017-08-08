function this = parameqapast(varargin)
%PARAMEQAPAST   Construct a PARAMEQAPAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.parameqapast;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
