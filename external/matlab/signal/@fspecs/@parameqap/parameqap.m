function this = parameqap(varargin)
%PARAMEQAP   Construct a PARAMEQAP object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.parameqap;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});


% [EOF]
