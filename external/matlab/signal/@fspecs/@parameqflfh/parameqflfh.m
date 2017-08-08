function this = parameqflfh(varargin)
%PARAMEQ   Construct a PARAMEQ object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.parameqflfh;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
