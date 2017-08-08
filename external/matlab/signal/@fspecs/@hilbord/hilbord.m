function this = hilbord(varargin)
%HILBORD   Construct a HILBORD object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.hilbord;

this.ResponseType = 'Hilbert Transformer with filter order';

this.FilterOrder = 30;

this.setspecs(varargin{:});

% [EOF]
