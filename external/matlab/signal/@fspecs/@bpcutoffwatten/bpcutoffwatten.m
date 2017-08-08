function this = bpcutoffwatten(varargin)
%BPCUTOFFWATTEN   Construct a BPCUTOFFWATTEN object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.bpcutoffwatten;

respstr = 'Bandpass with cutoff and attenuation';
fstart = 2;
fstop = 3;
nargsnoFs = 6;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});
