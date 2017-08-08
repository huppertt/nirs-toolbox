function this = hpstopapass(varargin)
%HPSTOPAPASS   Construct a HPSTOPAPASS object.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

this = fspecs.hpstopapass;

respstr = 'Highpass with passband ripple';
fstart = 2;
fstop = 2;
nargsnoFs = 4;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
