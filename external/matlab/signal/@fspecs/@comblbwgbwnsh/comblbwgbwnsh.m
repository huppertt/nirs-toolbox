function this = comblbwgbwnsh(varargin)
%COMBLBWGBWNSH   Construct a COMBLBWGBWNSH object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.comblbwgbwnsh;

set(this,'CombType','Notch');

set(this, 'ResponseType', 'Comb Filter');

this.setspecs(varargin{:});

% [EOF]
