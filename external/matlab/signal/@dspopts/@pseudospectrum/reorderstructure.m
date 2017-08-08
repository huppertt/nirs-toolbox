function s = reorderstructure(this, s) %#ok
%REORDERSTRUCTURE   

%   Author(s): W. Syed
%   Copyright 1988-2006 The MathWorks, Inc.

if (isprop(this, 'NFFT'))
    s = reorderstructure(s,'FreqPoints', 'NFFT','NormalizedFrequency','Fs','SpectrumRange', 'CenterDC');
elseif (isprop(this, 'FrequencyVector'))
    s = reorderstructure(s,'FreqPoints', 'FrequencyVector','NormalizedFrequency','Fs','SpectrumRange', 'CenterDC');
end


% [EOF]