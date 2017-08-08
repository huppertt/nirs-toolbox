function this = sbarbmagnphase(varargin)
%SBARBMAGNPHASE Construct a SBARBMAGNPHASE object.

%   Copyright 2005-2010 The MathWorks, Inc.

this = fspecs.sbarbmagnphase;

% Default response
set_defaultresponse(this)

respstr = 'Single-Band Arbitrary Magnitude and Phase';
fstart = 1;
fstop = 1;
nargsnoFs = 4;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});
% [EOF]
