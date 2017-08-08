function this = peakqap(varargin)
%PEAKQAP   Construct a PEAKQAP object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

this = fspecs.peakqap;

set(this, 'ResponseType', 'Peaking Filter');

this.setspecs(varargin{:});
% [EOF]
