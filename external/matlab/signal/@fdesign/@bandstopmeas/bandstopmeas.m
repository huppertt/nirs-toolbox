function this = bandstopmeas(hfilter, varargin)
%BANDSTOPMEAS   Construct a BANDSTOPMEAS object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

% Constructor an "empty" object.
this = fdesign.bandstopmeas;

% Parse the inputs for the fdesign object.
minfo = parseconstructorinputs(this, hfilter, varargin{:});

if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end

% Measure the bandstop filter remarkable frequencies.
this.Fpass1 = findfpass(this, reffilter(hfilter), minfo.Fpass1, minfo.Apass1, 'down', ...
    [0 minfo.Fcutoff1 minfo.Fstop1]);
this.F3dB1  = findfrequency(this, hfilter, 1/sqrt(2), 'down', 'first');
this.F6dB1  = findfrequency(this, hfilter, 1/2, 'down', 'first');
this.Fstop1 = findfstop(this, reffilter(hfilter), minfo.Fstop1, minfo.Astop, 'down');
this.Fstop2 = findfstop(this, reffilter(hfilter), minfo.Fstop2, minfo.Astop, 'up');
this.F6dB2  = findfrequency(this, hfilter, 1/2, 'up', 'last');
this.F3dB2  = findfrequency(this, hfilter, 1/sqrt(2), 'up', 'last');
this.Fpass2 = findfpass(this, reffilter(hfilter), minfo.Fpass2, minfo.Apass2, 'up', ...
    [max([minfo.Fstop2 minfo.Fcutoff2]) Fs/2]);

% Use the measured Fpass1, Fpass2, Fstop1 and Fstop2 when they are not
% specified to have a true measure of Apass1, Apass2 and Astop. See
% G425069.
if isempty(minfo.Fpass1), minfo.Fpass1 = this.Fpass1; end 
if isempty(minfo.Fpass2), minfo.Fpass2 = this.Fpass2; end 
if isempty(minfo.Fstop1), minfo.Fstop1 = this.Fstop1; end
if isempty(minfo.Fstop2), minfo.Fstop2 = this.Fstop2; end

% Measure ripples and attenuations.
this.Apass1 = measureripple(this, hfilter, 0, minfo.Fpass1, minfo.Apass1);
this.Astop  = measureattenuation(this, hfilter, minfo.Fstop1, minfo.Fstop2, minfo.Astop);
this.Apass2 = measureripple(this, hfilter, minfo.Fpass2, Fs/2, minfo.Apass2);

% [EOF]
