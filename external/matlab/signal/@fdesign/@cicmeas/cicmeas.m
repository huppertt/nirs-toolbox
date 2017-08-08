function this = cicmeas(hfilter, hfdesign)
%CICMEAS   Construct a CICMEAS object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct'));

this = fdesign.cicmeas;

% If the specification is not passed, use the stored specification.
if nargin < 2
    hfdesign = getfdesign(hfilter);
    if isempty(hfdesign)
        error(message('signal:fdesign:cicmeas:cicmeas:missingFDesign'));
    end
end

% Set CIC nulls
if isa(hfilter,'mfilt.cicdecim'),
    Fo = cicnulls(hfilter.DecimationFactor,hfilter.DifferentialDelay);
else
    Fo = cicnulls(hfilter.InterpolationFactor,hfilter.DifferentialDelay);
end

this.Fnulls = Fo;

% Sync up the sampling frequencies.
if ~hfdesign.NormalizedFrequency
    this.normalizefreq(false, get(hfdesign, 'fs'));
end

% Sync up Passband Frequency
this.Fpass = hfdesign.Fpass;

% Set stopband frequency
this.Fstop = this.Fnulls(1) - this.Fpass;

set(this, 'Specification', hfdesign);

% Find passband ripple and stopband attenuation
if this.NormalizedFrequency,
    H = freqz(hfilter,[0, this.Fpass, this.Fstop]*pi);
else
    H = freqz(hfilter,[0, this.Fpass, this.Fstop],this.Fs);
end

dcgain = 20*log10(abs(H(1)));

this.Apass = dcgain - 20*log10(abs(H(2)));
this.Astop = dcgain - 20*log10(abs(H(3)));

% [EOF]
