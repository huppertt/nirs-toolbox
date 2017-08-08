function s = reorderstructure(this, s) %#ok
%REORDERSTRUCTURE   

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

if (isprop(this, 'NFFT'))
    s = reorderstructure(s,'FreqPoints', 'NFFT','NormalizedFrequency','Fs','SpectrumType', 'CenterDC');
elseif (isprop(this, 'FrequencyVector'))
    s = reorderstructure(s,'FreqPoints', 'FrequencyVector','NormalizedFrequency','Fs','SpectrumType', 'CenterDC');
end


% [EOF]
