function this = parameqaudioshelfqa(varargin)
%PARAMEQ   Construct a PARAMEQAUDIOSHELFQA object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.parameqaudioshelfqa;

set(this, 'ResponseType', 'Parametric Equalizer');

this.setspecs(varargin{:});

% [EOF]
