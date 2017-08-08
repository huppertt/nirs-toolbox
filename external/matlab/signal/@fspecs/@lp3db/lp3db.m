function this = lp3db(varargin)
%LP3DB   Construct a LP3DB object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.lp3db;

set(this, 'ResponseType', 'Lowpass with 3-dB Frequency Point');

this.setspecs(varargin{:});

% [EOF]
