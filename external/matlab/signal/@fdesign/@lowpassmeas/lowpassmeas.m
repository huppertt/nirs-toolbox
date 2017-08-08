function this = lowpassmeas(hfilter, varargin)
%LOWPASSMEAS   Construct a LOWPASSMEAS object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

% Construct an "empty" object.
this = fdesign.lowpassmeas;

% Parse the inputs.
minfo = parseconstructorinputs(this, hfilter, varargin{:});

if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end
    
% Measure the lowpass filter remarkable frequencies.
this.Fpass = findfpass(this, reffilter(hfilter), minfo.Fpass, minfo.Apass, 'down');
this.F3dB  = findfrequency(this, hfilter, 1/sqrt(2), 'down', 'first');
this.F6dB  = findfrequency(this, hfilter, 1/2, 'down', 'first');
this.Fstop = findfstop(this, reffilter(hfilter), minfo.Fstop, minfo.Astop, 'down');

% Use the measured Fpass and Fstop when they are not specified to have a
% true measure of Apass and Astop. See G425069.
if isempty(minfo.Fpass), minfo.Fpass = this.Fpass; end 
if isempty(minfo.Fstop), minfo.Fstop = this.Fstop; end

% Measure ripples and attenuation.
this.Apass = measureripple(this, hfilter, 0, minfo.Fpass, minfo.Apass);
this.Astop = measureattenuation(this, hfilter, minfo.Fstop, Fs/2, minfo.Astop);

% [EOF]
