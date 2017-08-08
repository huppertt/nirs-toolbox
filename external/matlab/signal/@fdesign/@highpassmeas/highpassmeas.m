function this = highpassmeas(hfilter, varargin)
%HIGHPASSMEAS   Construct a HIGHPASSMEAS object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

% Construct and "empty" object.
this = fdesign.highpassmeas;

% Parse the inputs.
minfo = parseconstructorinputs(this, hfilter, varargin{:});

if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end

% Measure the highpass filter remarkable frequencies.
this.Fstop = findfstop(this, reffilter(hfilter), minfo.Fstop, minfo.Astop, 'up');
this.F6dB  = findfrequency(this, hfilter, 1/2, 'up', 'first');
this.F3dB  = findfrequency(this, hfilter, 1/sqrt(2), 'up', 'first');
this.Fpass = findfpass(this, reffilter(hfilter), minfo.Fpass, minfo.Apass, 'up');

% Use the measured Fpass and Fstop when they are not specified to have a
% true measure of Apass and Astop. See G425069.
if isempty(minfo.Fpass), minfo.Fpass = this.Fpass; end 
if isempty(minfo.Fstop), minfo.Fstop = this.Fstop; end

% Measure ripples and attenuation.
this.Astop = measureattenuation(this, hfilter, 0, minfo.Fstop, minfo.Astop);
this.Apass = measureripple(this, hfilter, minfo.Fpass, Fs/2, minfo.Apass);

% [EOF]
