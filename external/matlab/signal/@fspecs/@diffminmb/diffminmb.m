function this = diffminmb(varargin)
%DIFFMINMB   Construct a DIFFMINMB object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.diffminmb;

this.ResponseType = 'Minimum-order multi-band Differentiator';

this.setspecs(varargin{:});


% [EOF]
