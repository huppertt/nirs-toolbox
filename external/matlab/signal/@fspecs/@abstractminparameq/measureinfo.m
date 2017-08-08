function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

minfo.F0 = this.F0;
minfo.BW = this.BW;
if isprop(this,'bwpass') && ~(isprop(this,'gpass') && isprop(this,'gstop')),
    minfo.BWpass = this.BWpass;
else
    minfo.BWpass = [];
end
if isprop(this,'bwstop'),
    minfo.BWstop = this.BWstop;
else
    minfo.BWstop = [];
end
if this.NormalizedFrequency,
    [Flow,Fhigh] = parameqbandedge(minfo.F0*pi,minfo.BW*pi,0);
    minfo.Flow   = Flow/pi;
    minfo.Fhigh  = Fhigh/pi;
else
    Fs = this.Fs;
    [Flow,Fhigh] = parameqbandedge(pi*minfo.F0/(Fs/2),pi*minfo.BW/(Fs/2),0);
    minfo.Flow   = Flow*Fs/(2*pi);
    minfo.Fhigh  = Fhigh*Fs/(2*pi);
end
minfo.Gref   = [];
minfo.G0     = [];
minfo.GBW    = this.GBW;
if isprop(this,'gpass') && isprop(this,'gstop'),
    minfo.Gstop = this.Gstop;
    minfo.Gpass = this.Gpass;
else
    % Determine through measurements
    minfo.Gstop = [];
    minfo.Gpass = [];
end


% [EOF]
