function this = hpcutoffwap(varargin)
%HPCUTOFFWAP   Construct a HPCUTOFFWAP object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.hpcutoffwap;

respstr = 'Highpass with cutoff and passband ripple';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
