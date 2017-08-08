function ylbl = getylabel(this)
%GETYLABEL Get the ylabel.

%   Copyright 2007 The MathWorks, Inc.

if this.NormalizedFrequency
    ylbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBradsample'));
else
    ylbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBHz'));
end


% [EOF]
