function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Copyright 2008 The MathWorks, Inc.

minfo.F0 = this.F0;

minfo.BWpass = [];
minfo.BWstop = [];
if this.NormalizedFrequency,
    minfo.BW = 2*atan(sin(pi*this.F0)/(2*this.Qa))/pi;
    [Flow,Fhigh] = parameqbandedge(minfo.F0*pi,minfo.BW*pi,0);
    minfo.Flow   = Flow/pi;
    minfo.Fhigh  = Fhigh/pi;
else
    Fs = this.Fs;
    w0 = pi*minfo.F0/(Fs/2);
    minfo.BW = (2*atan(sin(w0)/(2*this.Qa))/pi)*(Fs/2);
    [Flow,Fhigh] = parameqbandedge(pi*minfo.F0/(Fs/2),pi*minfo.BW/(Fs/2),0);
    minfo.Flow   = Flow*Fs/(2*pi);
    minfo.Fhigh  = Fhigh*Fs/(2*pi);
end
minfo.Gref   = [];
minfo.G0     = [];
minfo.GBW = 10*log10(sqrt(10^(this.G0/10)* 10^(this.Gref/10)));
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
