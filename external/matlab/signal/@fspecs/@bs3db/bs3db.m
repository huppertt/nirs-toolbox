function this = bs3db(varargin)
%BS3DB   Construct a BS3DB object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.bs3db;

set(this, 'ResponseType', 'Bandstop with 3-dB Frequency Point');

this.setspecs(varargin{:});

% [EOF]
