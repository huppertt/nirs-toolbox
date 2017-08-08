function this = peakbwap(varargin)
%PEAKBWAP   Construct a PEAKBWAP object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.peakbwap;

set(this, 'ResponseType', 'Peaking Filter');

this.setspecs(varargin{:});

% [EOF]
