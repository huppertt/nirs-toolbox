function thisnormalizefreq(this,oldFs,newFsFlag)
%THISNORMALIZEFREQ   Normalize/un-normalize the frequency of the data object.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

% Scale the data by the appropriate value, either Fs, oldFs, or 2pi.
Fs = getprivfs(this);
if this.NormalizedFrequency,  
    scaleFactor = oldFs/(2*pi);
else
    if newFsFlag,  % Catch case of repeated calls with new Fs.
        scaleFactor = oldFs/Fs;
    else
        scaleFactor = (2*pi)/Fs;
    end
end
this.Data = this.Data*scaleFactor;

% [EOF]
