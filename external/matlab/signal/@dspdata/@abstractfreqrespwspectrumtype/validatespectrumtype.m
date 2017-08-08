function validatespectrumtype(this,spectrumType)
%VALIDATESPECTRUMTYPE   Validate SpectrumType property value.
%
% This error checking should be done in the object's set method, but for
% enum datatypes UDD first checks the list before calling the set method.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

validStrs = {'onesided','twosided'};
if ~ischar(spectrumType) | ~any(strcmpi(spectrumType,validStrs)),
    error(message('signal:dspdata:abstractfreqrespwspectrumtype:validatespectrumtype:invalidSpectrumType', 'SpectrumType', validStrs{ 1 }, validStrs{ 2 }));
end


% [EOF]
