function this = parameqaudioshelfs(varargin)
%PARAMEQ   Construct a PARAMEQAUDIOSHELFS object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.parameqaudioshelfs;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
