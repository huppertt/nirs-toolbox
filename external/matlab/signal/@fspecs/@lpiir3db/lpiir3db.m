function this = lpiir3db(varargin)
%LP3DB   Construct a LPIIR3DB object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.lpiir3db;

constructor(this,varargin{:});

this.ResponseType = 'Lowpass with 3-dB Frequency Point';

% [EOF]