function this = combq(varargin)
%COMBBW   Construct a COMBQ object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.combq;

set(this, 'ResponseType', 'Comb Filter');

this.setspecs(varargin{:});

% [EOF]
