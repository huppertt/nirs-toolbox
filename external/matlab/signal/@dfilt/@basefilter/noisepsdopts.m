function opts = noisepsdopts(this)
%NOISEPSDOPTS   Options for noisepsd method.

%   Author(s): R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.

% Construct default opts object
opts = dspopts.spectrum;

% Make sure NFFT is not too large, this method takes a while
opts.NFFT = 512;

if isreal(this),
    opts.SpectrumType = 'Onesided';
else
    opts.SpectrumType = 'Twosided';
end

% Don't center DC by default
opts.CenterDC = false;

% [EOF]
