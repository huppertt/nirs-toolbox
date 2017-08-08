function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Copyright 2008-2011 The MathWorks, Inc.

minfo.Fpass = 2*this.BT/this.SamplesPerSymbol;
% Note that:
% Fpass is in normalized frequency
% Fpass = B/(Fs/2);  Fs is the sampling frequency
% Fpass = B/(Fs/2) * T/(sps/Fs);
% Fpass = 2*B*T/sps;
minfo.Fstop = [];
minfo.Apass = [];
minfo.Astop = [];

if ~this.NormalizedFrequency
  minfo.Fpass = minfo.Fpass*this.Fs/2;
  minfo.Fstop = minfo.Fstop*this.Fs/2;
end

% [EOF]
