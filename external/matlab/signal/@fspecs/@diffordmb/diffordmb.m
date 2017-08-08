function this = diffordmb(varargin)
%DIFFORDMB   Construct a DIFFORDMB object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.diffordmb;

this.ResponseType = 'Multi-band Differentiator with filter order';

% Defaults
this.FilterOrder = 30;
this.Fpass = .7;
this.Fstop = .9;  

this.setspecs(varargin{:});

% [EOF]
