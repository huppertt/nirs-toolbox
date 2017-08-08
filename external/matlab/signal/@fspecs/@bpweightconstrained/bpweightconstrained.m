function this = bpweightconstrained(varargin)
%BPWEIGHTCONSTRAINED Construct a BPWEIGHTCONSTRAINED object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.bpweightconstrained;

respstr = 'Bandpass';
fstart = 2;
fstop = 5;
nargsnoFs = 8;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
