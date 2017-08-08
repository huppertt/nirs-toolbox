function this = sbarbmagmin(varargin)
%SBARBMAGMIN Construct a SBARBMAGMIN object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.sbarbmagmin;

respstr = 'Single-Band Arbitrary Magnitude';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});
% [EOF]
