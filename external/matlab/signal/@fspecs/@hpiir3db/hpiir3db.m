function this = hpiir3db(varargin)
%LP3DB Construct a HPIIR3DB object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.hpiir3db;

constructor(this,varargin{:});

this.ResponseType = 'Highpass with 3-dB Frequency Point';

% [EOF]