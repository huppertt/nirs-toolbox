function s = aswfs_getdesignpanelstate(this)
%ASWFS_GETDESIGNPANELSTATE   

%   Author(s): J. Schickler
%   Copyright 2004-2005 The MathWorks, Inc.

if this.NormalizedFrequency
    s.Components{1}.freqUnits = 'Normalized (0 to 1)';
else
    s.Components{1}.freqUnits = 'Hz';
    s.Components{1}.Fs        = sprintf('%d', this.Fs);
end

% [EOF]
