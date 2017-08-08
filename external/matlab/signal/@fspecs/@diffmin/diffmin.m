function this = diffmin(varargin)
%DIFFMIN   Construct a DIFFMIN object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.diffmin;

this.ResponseType = 'Minimum-order Differentiator';

this.setspecs(varargin{:});

% [EOF]
