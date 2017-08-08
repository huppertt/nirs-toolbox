function this = peakqast(varargin)
%PEAKQAST   Construct a PEAKQAST object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.peakqast;

set(this, 'ResponseType', 'Peaking Filter');

this.setspecs(varargin{:});

% [EOF]
