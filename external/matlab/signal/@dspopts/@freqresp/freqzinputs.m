function c = freqzinputs(this)
%FREQZINPUTS   Return a cell with the inputs for FREQZ, PHASEZ, etc.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

switch lower(this.FrequencySpecification)
    case 'nfft'
        c = {this.NFFT, this.SpectrumRange};
    case 'frequencyvector'
        c = {this.FrequencyVector};
end

if ~this.NormalizedFrequency
    c = {c{:}, this.Fs};
end

% [EOF]
