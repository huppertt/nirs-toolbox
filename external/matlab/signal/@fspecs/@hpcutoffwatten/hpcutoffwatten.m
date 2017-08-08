function this = hpcutoffwatten(varargin)
%HPCUTOFFWATTEN   Construct a HPCUTOFFWATTEN object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.hpcutoffwatten;

respstr = 'Highpass with cutoff and attenuation';
fstart = 2;
fstop = 2;
nargsnoFs = 4;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
