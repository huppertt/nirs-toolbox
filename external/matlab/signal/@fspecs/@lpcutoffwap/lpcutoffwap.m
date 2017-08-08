function this = lpcutoffwap(varargin)
%LPCUTOFFWAP   Construct a LPCUTOFFWAP object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.lpcutoffwap;

respstr = 'Lowpass with cutoff and passband ripple';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
