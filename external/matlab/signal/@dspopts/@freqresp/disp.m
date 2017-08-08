function disp(this)
%DISP   Display this object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

p = {'NormalizedFrequency'};

if ~this.NormalizedFrequency
    p = {p{:}, 'Fs'};
end

p = {p{:}, 'FrequencySpecification'};

if strcmpi(this.FrequencySpecification, 'NFFT')
    p = {p{:}, 'NFFT', 'SpectrumRange', 'CenterDC'};
else
    p = {p{:}, 'FrequencyVector'};
end

siguddutils('dispstr', this, p);

% [EOF]
