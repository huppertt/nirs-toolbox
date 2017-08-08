function this = combbw(varargin)
%COMBBW   Construct a COMBBW object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.combbw;

set(this, 'ResponseType', 'Comb Filter');

this.setspecs(varargin{:});

% [EOF]
