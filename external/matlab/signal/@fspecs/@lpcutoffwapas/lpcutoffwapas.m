function this = lpcutoffwapas(varargin)
%LPCUTOFFWAPAS   Construct a LPCUTOFFWAPAS object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.lpcutoffwapas;

respstr = 'Lowpass with cutoff, passband ripple and stopband attenuation';
fstart = 2;
fstop = 2;
nargsnoFs = 4;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
