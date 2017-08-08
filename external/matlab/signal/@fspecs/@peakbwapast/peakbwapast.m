function this = peakbwapast(varargin)
%PEAKBWAPAST   Construct a PEAKBWAPAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.peakbwapast;

set(this, 'ResponseType', 'Peaking Filter');

this.setspecs(varargin{:});

% [EOF]
