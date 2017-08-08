function this = bpcutoffwapas(varargin)
%BPCUTOFFWAPAS   Construct a BPCUTOFFWAPAS object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.bpcutoffwapas;

respstr = 'Bandpass with cutoff, passband ripple and stopband attenuation';
fstart = 2;
fstop = 2;
nargsnoFs = 4;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
