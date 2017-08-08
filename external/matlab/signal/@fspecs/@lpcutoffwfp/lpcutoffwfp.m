function this = lpcutoffwfp(varargin)
%LPCUTOFFWFP   Construct a LPCUTOFFWFP object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.lpcutoffwfp;

respstr = 'Lowpass with cutoff and passband frequency';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
