function this = isinclpmeas(hfilter, ~, varargin)
%ISINCLPMEAS Construct an ISINCLPMEAS object.

%   Copyright 2005-2011 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

this = fdesign.isinclpmeas;

% Parse the inputs.
minfo = parseconstructorinputs(this, hfilter, varargin{:});
    
if this.NormalizedFrequency, Fs = 2;
else                        Fs = this.Fs; end

idealfcn = {@defineideal, minfo.FrequencyFactor, ...
     minfo.Power, minfo.CICRateChangeFactor, ...
    [minfo.Fpass/(Fs/2) minfo.Fcutoff/(Fs/2) minfo.Fstop/(Fs/2)]};

% Normalize numerator according to interpolation factor  
hfilter1 = copy(hfilter);
rcf = getratechangefactors(hfilter);
hfilter1.Numerator = hfilter1.Numerator/rcf(1);
  
% Measure the lowpass filter remarkable frequencies.
this.Fpass = findfpass(this, reffilter(hfilter1), minfo.Fpass, minfo.Apass, 'down', ...
    [0 Fs/2], idealfcn);
this.F3dB  = findfrequency(this, hfilter, 1/sqrt(2), 'down', 'first');
this.F6dB  = findfrequency(this, hfilter, 1/2, 'down', 'first');
this.Fstop = findfstop(this, reffilter(hfilter), minfo.Fstop, minfo.Astop, 'down');

% Use the measured Fpass and Fstop when they are not specified to have a
% true measure of Apass and Astop. See G425069.
if isempty(minfo.Fpass), minfo.Fpass = this.Fpass; end 
if isempty(minfo.Fstop), minfo.Fstop = this.Fstop; end

% Measure ripples and attenuation.
this.Apass = measureripple(this, hfilter1, 0, minfo.Fpass, minfo.Apass, idealfcn);
this.Astop = measureattenuation(this, hfilter, minfo.Fstop, Fs/2, minfo.Astop);

% -------------------------------------------------------------------------
function H = defineideal(F, FreqFactor, Power, CICRCF, Fpass)

Fpass = Fpass(1);

if CICRCF == 1
  H = 1./sinc(F*FreqFactor).^Power;
else
  H = computeInvDirichletSinc(F, FreqFactor, Power, CICRCF);
end

% Make sure that we don't apply the inverse sinc into the stopband.
indx = find(F > Fpass, 1, 'first');

if ~isempty(indx)
    H(indx+1:end) = H(indx);
end

% ---------------------------------------------------
function y = computeInvDirichletSinc(f, FreqFactor, Power, CICRCF)

  % CIC differential delay is equal to the frequency factor times 2
  diffDelay = FreqFactor*2;
  x = f*pi;
  y = ones(size(x));  
  i = find(sin(x/(2*CICRCF)) & sin(diffDelay*x/2));  

  y(i) = (CICRCF*diffDelay * sin(x(i)/(2*CICRCF))./sin(diffDelay*x(i)/2)).^Power;  

% [EOF]
