function this = hpcutoffwfp(varargin)
%HPCUTOFFWFP   Construct a HPCUTOFFWFP object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.hpcutoffwfp;

respstr = 'Highpass with cutoff and passband frequency';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
