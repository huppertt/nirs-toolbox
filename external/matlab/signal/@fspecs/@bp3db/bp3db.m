function this = bp3db(varargin)
%BP3DB   Construct a BP3DB object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.bp3db;

set(this, 'ResponseType', 'Bandpass with 3-dB Frequency Point');

this.setspecs(varargin{:});

% [EOF]
