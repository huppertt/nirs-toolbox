function this = lpstopapass(varargin)
%LPSTOPAPASS   Construct a LPSTOPAPASS object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

this = fspecs.lpstopapass;

respstr = 'Lowpass with passband ripple';
fstart = 1;
fstop = 2;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
