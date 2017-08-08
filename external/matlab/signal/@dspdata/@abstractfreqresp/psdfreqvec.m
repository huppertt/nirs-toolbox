function W = psdfreqvec(this,hopts) %#ok<INUSL>
%PSDFREQVEC  Returns the frequency vector with appropriate frequency range.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

% Determine length of FFT to calculate the spectrum over the whole Nyquist
% interval.
lenX = hopts.NFFT;
wholenfft = lenX;
range = 'whole';

if ishalfnyqinterval(hopts);
    range = 'half';
    wholenfft = (lenX-1)*2; % Spec: always assume whole was EVEN; required by psdfreqvec.m
end

% psdfreqvec requires the whole NFFT to calculate 'half' Nyquist range.
if hopts.NormalizedFrequency,
    Fs = []; % Used by psdfreqvec to calculate rad/sample.
else
    Fs = hopts.Fs;
end

% Special case the singleton/DC case.
if wholenfft == 0
    W = 0;
else
    W = psdfreqvec('npts',wholenfft,'Fs',Fs,'CenterDC',hopts.CenterDC,'Range',range);
end

% [EOF]
