function this = bpcutoffwap(varargin)
%BPCUTOFFWAP   Construct a BPCUTOFFWAP object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.bpcutoffwap;

respstr = 'Bandpass with cutoff and passband ripple';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
