function this = cicinterpspecs
%CICINTERPSPECS   Construct a CICINTERPSPECS object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.cicinterpspecs;

this.ResponseType = 'CIC Interpolator';

% Set the default normalized passband frequency to a small value which will
% result in a design of N = 2
this.Fpass = .05;

% [EOF]
