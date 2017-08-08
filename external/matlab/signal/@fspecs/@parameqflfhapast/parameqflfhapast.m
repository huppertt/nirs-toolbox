function this = parameqflfhapast(varargin)
%PARAMEQFLFHAPAST   Construct a PARAMEQFLFHAPAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.parameqflfhapast;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
