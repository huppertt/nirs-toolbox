function Fs = setfs(h,Fs)
%SETFS   

%   Copyright 1988-2011 The MathWorks, Inc.

if Fs == 0,
    error(message('signal:fspecs:abstractspecwithfs:setfs:invalidSpec'));
end

if Fs == Inf,
    error(message('signal:fspecs:abstractspecwithfs:setfs:invalidSpecInf'));
end

if h.NormalizedFrequency,
    error(message('signal:fspecs:abstractspecwithfs:setfs:readonly', 'Fs', 'NormalizedFrequency', 'normalizefreq(h,false,Fs)'));        
end

h.privFs = Fs;

% Make Fs empty to not duplicate storage
Fs = [];

% [EOF]
