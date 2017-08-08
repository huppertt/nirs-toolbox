function s = reorderstructure(this,s) %#ok
%REORDERSTRUCTURE   

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

s = reorderstructure(s,'NFFT','NormalizedFrequency','Fs','SpectrumRange', 'CenterDC');

% [EOF]
