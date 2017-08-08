function this = bscutoffwbws(varargin)
%BSCUTOFFWBWS   Construct a BSCUTOFFWBWS object.

%   Author(s): V. Pellissier
%   Copyright 2004 The MathWorks, Inc.

this = fspecs.bscutoffwbws;

respstr = 'Bandstop with cutoff and stopband width';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
