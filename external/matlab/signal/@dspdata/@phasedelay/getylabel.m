function ylbl = getylabel(this)
%GETYLABEL   Get the ylabel.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

if this.NormalizedFrequency
    ylbl = getString(message('signal:dspdata:dspdata:PhaseDelaySamples'));
else
    ylbl = getString(message('signal:dspdata:dspdata:PhaseDelayRadiansHz'));
end

% [EOF]
