function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.
if this.NormalizedFrequency,
    [F0,BW] = parameqbandedge(this.Flow*pi,this.Fhigh*pi,1);
    minfo.F0 = F0/pi;
    minfo.BW = BW/pi;
else
    Fs = this.Fs;
    [F0,BW] = parameqbandedge(pi*this.Flow/(Fs/2),pi*this.Fhigh/(Fs/2),1);
    minfo.F0 = F0*Fs/(2*pi);
    minfo.BW = BW*Fs/(2*pi);
end

minfo.BWpass = [];
minfo.BWstop = [];
minfo.Flow   = this.Flow;
minfo.Fhigh  = this.Fhigh;
minfo.Gref   = [];
minfo.G0     = [];
minfo.GBW    = this.GBW;
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
