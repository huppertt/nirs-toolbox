function this = diffordmbap(varargin)
%DIFFORDMBAP Construct a DIFFORDMBAP object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.diffordmbap;

this.ResponseType = 'Multi-band Differentiator with filter order';

% Defaults
this.FilterOrder = 30;
this.Fpass = .7;
this.Fstop = .9;  
this.Apass = 1;

this.setspecs(varargin{:});

% [EOF]
