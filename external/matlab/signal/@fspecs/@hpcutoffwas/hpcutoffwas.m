function this = hpcutoffwas(varargin)
%HPCUTOFFWAS   Construct a HPCUTOFFWAS object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.hpcutoffwas;

respstr = 'Highpass with cutoff and stopband ripple';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
