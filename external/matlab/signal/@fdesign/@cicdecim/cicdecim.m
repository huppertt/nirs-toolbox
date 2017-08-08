function this = cicdecim(varargin)
%CICDECIM   Construct a CICDECIM object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

this = fdesign.cicdecim;

set(this, 'Response', 'CIC Decimator');

this.setspecs(varargin{:});

% [EOF]
