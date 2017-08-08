function this = bscutoffwas(varargin)
%BSCUTOFFWAS   Construct a BSCUTOFFWAS object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.bscutoffwas;

respstr = 'Bandstop with cutoff and stopband attenuation';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
