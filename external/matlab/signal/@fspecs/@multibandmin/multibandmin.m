function this = multibandmin(varargin)
%MULTIBANDMIN Construct a MULTIBANDMIN object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.multibandmin;

respstr = 'Multi-Band Arbitrary Magnitude';
fstart = 1;
fstop = 1;
nargsnoFs = 2;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
