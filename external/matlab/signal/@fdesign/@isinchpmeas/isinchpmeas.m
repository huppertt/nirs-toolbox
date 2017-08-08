function this = isinchpmeas(hfilter, ~, varargin)
%ISINCHPMEAS Construct an ISINCHPMEAS object.

%   Copyright 2011 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

this = fdesign.isinchpmeas;

% Parse the inputs.
minfo = parseconstructorinputs(this, hfilter, varargin{:});
    
if this.NormalizedFrequency, Fs = 2;
else                        Fs = this.Fs; end

idealfcn = {@defineideal, minfo.FrequencyFactor, minfo.Power, ...
    [minfo.Fpass/(Fs/2) minfo.Fcutoff/(Fs/2) minfo.Fstop/(Fs/2)]};

% Normalize numerator according to interpolation factor  
hfilter1 = copy(hfilter);
rcf = getratechangefactors(hfilter);
hfilter1.Numerator = hfilter1.Numerator/rcf(1);

% Measure the highpass filter remarkable frequencies.
this.Fstop = findfstop(this, reffilter(hfilter), minfo.Fstop, minfo.Astop, 'up');
this.F6dB  = findfrequency(this, hfilter, 1/2, 'up', 'first');
this.F3dB  = findfrequency(this, hfilter, 1/sqrt(2), 'up', 'first');
this.Fpass = findfpass(this, reffilter(hfilter1), minfo.Fpass, minfo.Apass, 'up',[0 Fs/2], idealfcn);
  
% Use the measured Fpass and Fstop when they are not specified to have a
% true measure of Apass and Astop.
if isempty(minfo.Fpass), minfo.Fpass = this.Fpass; end 
if isempty(minfo.Fstop), minfo.Fstop = this.Fstop; end

% Measure ripples and attenuation.
this.Astop = measureattenuation(this, hfilter, 0, minfo.Fstop, minfo.Astop);
this.Apass = measureripple(this, hfilter1, minfo.Fpass, Fs/2, minfo.Apass, idealfcn);

% -------------------------------------------------------------------------
function H = defineideal(F, FreqFactor, Power, Fpass)

Fpass = Fpass(1);

H = 1./sinc((1-F)*FreqFactor).^Power;

% Make sure that we don't apply the inverse sinc in the stopband.
indx = find(F < Fpass, 1, 'first');

if ~isempty(indx)
    H(indx+1:end) = H(indx);
end

% [EOF]
