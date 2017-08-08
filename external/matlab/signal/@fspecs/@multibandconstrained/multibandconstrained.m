function this = multibandconstrained(varargin)
%MULTIBANDCONSTRAINED Construct a MULTIBANDCONSTRAINED object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.multibandconstrained;

respstr = 'Multi-Band Arbitrary Magnitude';
fstart = 1;
fstop = 1;
nargsnoFs = 2;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
