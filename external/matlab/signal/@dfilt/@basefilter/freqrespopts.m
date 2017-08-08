function opts = freqrespopts(this)
%FREQRESPOPTS  Create an options object for frequency response estimate.
%   OPTS = FREQRESPOPTS(Hd) uses the current settings in the filter Hd to
%   create an options object OPTS that contains options for frequency
%   response estimation.  The OPTS object can be passed in as an argument
%   to the FREQRESPEST method.
%  
%   See also DFILT/FREQRESPEST, DFILT/NOISEPSD, DFILT/NOISEPSDOPTS,
%   DFILT/SCALE,  DFILT/NORM.  

%   Author(s): R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.

% Construct default opts object
opts = dspopts.pseudospectrum;

% Make sure NFFT is not too large, this method takes a while
opts.NFFT = 512;

if isreal(this),
    opts.SpectrumRange = 'Half';
else
    opts.SpectrumRange = 'Whole';
end

% Don't center DC by default
opts.CenterDC = false;


% [EOF]
