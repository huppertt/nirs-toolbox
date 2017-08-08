function this = lpstopfpass(varargin)
%LPSTOPFPASS   Construct a LPSTOPFPASS object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

this = fspecs.lpstopfpass;

respstr = 'Lowpass with passband frequency';
fstart = 1;
fstop = 2;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
