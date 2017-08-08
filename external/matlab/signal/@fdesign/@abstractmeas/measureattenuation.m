function atten = measureattenuation(this, hfilter, Fstart, Fend, Astop)
%MEASUREATTENUATION   Return the attenuation in the stopband.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(4,5,nargin,'struct'));

if nargin < 5
    Astop = [];
end

hfdesign = get(this, 'Specification');

if isempty(Fstart) || isempty(Fend)
    atten = Astop;
else

    N = 2^12;

    % calculate max of response from w_lo to w_hi
    if this.NormalizedFrequency
        Fs = 2;
    else
        Fs = this.Fs;
    end

    % Calculate the frequency response in the stopband.
    h = abs(freqz(hfilter, linspace(Fstart, Fend, N), Fs));

    % The attenuation is defined as the distance between the nominal gain
    % and the maximum rippple in the stopband.
    ngain = nominalgain(hfilter);
    if isempty(ngain)
        ngain = 1;
    end
    atten = db(ngain)-db(max(h));
end

% [EOF]
