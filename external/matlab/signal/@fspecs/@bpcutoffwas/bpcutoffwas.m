function this = bpcutoffwas(varargin)
%BPCUTOFFWAS   Construct a BPCUTOFFWAS object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.bpcutoffwas;

respstr = 'Bandpass with cutoff and stopband attenuation';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
