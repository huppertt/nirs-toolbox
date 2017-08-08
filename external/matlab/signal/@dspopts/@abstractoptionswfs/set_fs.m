function Fs = set_fs(h,Fs)
%SETFS   

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
   
h.privFs = Fs;

% Unset NormalizedFrequency
h.NormalizedFrequency = false;

% Make Fs empty to not duplicate storage
Fs = [];

% [EOF]
