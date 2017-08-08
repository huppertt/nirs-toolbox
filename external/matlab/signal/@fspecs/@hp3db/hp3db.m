function this = hp3db(varargin)
%HP3DB   Construct a HP3DB object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.hp3db;

set(this, 'ResponseType', 'Highpass with 3-dB Frequency Point');

this.setspecs(varargin{:});

% [EOF]
