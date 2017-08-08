function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Copyright 2008 The MathWorks, Inc.
if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end

minfo.F0 = this.F0/(Fs/2);
if minfo.F0
    minfo.BW = 1-this.Fc/(Fs/2);
else
    minfo.BW = this.Fc/(Fs/2);
end
minfo.BWpass = [];
minfo.BWstop = [];

[Flow,Fhigh] = parameqbandedge(minfo.F0*pi,minfo.BW*pi,0);
minfo.Flow   = (Flow/pi)*(Fs/2);
minfo.Fhigh  = (Fhigh/pi)*(Fs/2);
minfo.BW = minfo.BW*(Fs/2);
minfo.F0 = minfo.F0*(Fs/2);

minfo.Gref   = [];
minfo.G0     = [];
minfo.GBW    = this.G0/2;
if isprop(this,'gpass'),
    minfo.Gpass = this.Gpass;
else
    minfo.Gpass = [];
end
if isprop(this,'gstop'),
    minfo.Gstop = this.Gstop;
else
    minfo.Gstop = [];
end

% [EOF]
