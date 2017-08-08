function ylbl = getylabel(this)
%GETYLABEL Get the ylabel.

%   Copyright 2007 The MathWorks, Inc.

if this.plotindb             
    ylbl = getString(message('signal:dspdata:dspdata:MagnitudedB'));
else
    ylbl = getString(message('signal:dspdata:dspdata:Magnitude'));
end


% [EOF]
