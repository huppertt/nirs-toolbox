function this = peakbw(varargin)
%PEAKBW   Construct a PEAKBW object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.peakbw;

set(this, 'ResponseType', 'Peaking Filter');

this.setspecs(varargin{:});

% [EOF]
