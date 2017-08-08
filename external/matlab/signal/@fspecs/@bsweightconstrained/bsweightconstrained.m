function this = bsweightconstrained(varargin)
%BSWEIGHTCONSTRAINED Construct a BSWEIGHTCONSTRAINED object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.bsweightconstrained;

respstr = 'Bandstop';
fstart = 2;
fstop = 5;
nargsnoFs = 8;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
