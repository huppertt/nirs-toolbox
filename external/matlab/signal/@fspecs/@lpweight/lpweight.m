function this = lpweight(varargin)
%LPWEIGHT   Construct a LPWEIGHT object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

this = fspecs.lpweight;

respstr = 'Lowpass';
fstart = 2;
fstop = 3;
nargsnoFs = 5;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
