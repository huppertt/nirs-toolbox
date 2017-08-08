function this = hpisinccutoffwatten(varargin)
%HPISINCCUTOFFWATTEN Construct a HPISINCCUTOFFWATTEN object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.hpisinccutoffwatten;

fstart = 2;
fstop = 2;
nargsnoFs = 6;
fsconstructor(this,'Inverse-sinc highpass',fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
