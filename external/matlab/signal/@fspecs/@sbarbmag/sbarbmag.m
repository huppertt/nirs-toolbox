function this = sbarbmag(varargin)
%SBARBMAG   Construct a SBARBMAG object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.sbarbmag;

respstr = 'Single-Band Arbitrary Magnitude';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});
% [EOF]
