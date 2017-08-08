function this = bandpassmeas(hfilter, varargin)
%BANDPASSMEAS   Construct a BANDPASSMEAS object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(1, inf, nargin,'struct'));

% Create the "empty" object.
this = fdesign.bandpassmeas;

% Parse the inputs to get the specification and the measurements list.
minfo = parseconstructorinputs(this, hfilter, varargin{:});

if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end

% Measure the bandpass filter remarkable frequencies.
this.Fstop1 = findfstop(this, reffilter(hfilter), minfo.Fstop1, minfo.Astop1, 'up', ...
    [0 minfo.Fpass1 minfo.Fcutoff1]);
this.F6dB1  = findfrequency(this, hfilter, 1/2, 'up', 'first');
this.F3dB1  = findfrequency(this, hfilter, 1/sqrt(2), 'up', 'first');
this.Fpass1 = findfpass(this, reffilter(hfilter), minfo.Fpass1, minfo.Apass, 'up');
this.Fpass2 = findfpass(this, reffilter(hfilter), minfo.Fpass2, minfo.Apass, 'down');
this.F3dB2  = findfrequency(this, hfilter, 1/sqrt(2), 'down', 'last');
this.F6dB2  = findfrequency(this, hfilter, 1/2, 'down', 'last');
this.Fstop2 = findfstop(this, reffilter(hfilter), minfo.Fstop2, minfo.Astop2, 'down', ...
    [max([minfo.Fpass2, minfo.Fcutoff2]) Fs/2]);

% Use the measured Fpass1, Fpass2, Fstop1 and Fstop2 when they are not
% specified to have a true measure of Apass, Astop1 and Astop2. See
% G425069.
if isempty(minfo.Fpass1), minfo.Fpass1 = this.Fpass1; end 
if isempty(minfo.Fpass2), minfo.Fpass2 = this.Fpass2; end 
if isempty(minfo.Fstop1), minfo.Fstop1 = this.Fstop1; end
if isempty(minfo.Fstop2), minfo.Fstop2 = this.Fstop2; end

% Measure ripples and attenuations.
this.Astop1 = measureattenuation(this, hfilter, 0, minfo.Fstop1, minfo.Astop1);
this.Apass  = measureripple(this, hfilter, minfo.Fpass1, minfo.Fpass2, minfo.Apass);
this.Astop2 = measureattenuation(this, hfilter, minfo.Fstop2, Fs/2, minfo.Astop2);

% [EOF]
