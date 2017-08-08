function this = parameqaudioqa(varargin)
%PARAMEQ   Construct a PARAMEQAUDIOQA object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.parameqaudioqa;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
