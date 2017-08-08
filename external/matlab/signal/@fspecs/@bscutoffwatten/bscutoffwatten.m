function this = bscutoffwatten(varargin)
%BSCUTOFFWATTEN   Construct a BSCUTOFFWATTEN object.

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.bscutoffwatten;

respstr = 'Bandstop with cutoff and attenuation';
fstart = 2;
fstop = 3;
nargsnoFs = 6;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});
