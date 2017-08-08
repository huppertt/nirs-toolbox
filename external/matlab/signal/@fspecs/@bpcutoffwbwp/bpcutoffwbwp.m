function this = bpcutoffwbwp(varargin)
%BPCUTOFFWBWP   Construct a BPCUTOFFWBWP object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.bpcutoffwbwp;

respstr = 'Bandpass with cutoff and passband width';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
