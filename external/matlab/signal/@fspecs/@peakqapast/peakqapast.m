function this = peakqapast(varargin)
%PEAKQAPAST   Construct a PEAKQAPAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.peakqapast;

set(this, 'ResponseType', 'Peaking Filter');

this.setspecs(varargin{:});

% [EOF]
