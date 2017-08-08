function this = cicinterp(varargin)
%CICINTERP   Construct a CICINTERP object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

this = fdesign.cicinterp;

set(this, 'Response', 'CIC Interpolator');

this.setspecs(varargin{:});

% [EOF]
