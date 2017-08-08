function this = lpcutoffwfs(varargin)
%LPCUTOFFWFS   Construct a LPCUTOFFWFS object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.lpcutoffwfs;

respstr = 'Lowpass with cutoff and stopband frequency';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% % [EOF]
