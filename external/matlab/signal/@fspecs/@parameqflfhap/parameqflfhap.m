function this = parameqflfhap(varargin)
%PARAMEQFLFHAP   Construct a PARAMEQFLFHAP object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.parameqflfhap;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
