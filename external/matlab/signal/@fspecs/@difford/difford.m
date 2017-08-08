function this = difford(varargin)
%DIFFORD   Construct a DIFFORD object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.difford;

this.ResponseType = 'Differentiator with filter order';

% Since this specification type can only be used to design type IV
% differentiators, set the default to an odd filter order.
this.FilterOrder = 31;

this.setspecs(varargin{:});


% [EOF]
