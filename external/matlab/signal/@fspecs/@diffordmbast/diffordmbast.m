function this = diffordmbast(varargin)
%DIFFORDMBAST Construct a DIFFORDMBAST object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.diffordmbast;

this.ResponseType = 'Multi-band Differentiator with filter order';

% Defaults
this.FilterOrder = 30;
this.Fpass = .7;
this.Fstop = .9;  
this.Astop = 60;

this.setspecs(varargin{:});

% [EOF]
