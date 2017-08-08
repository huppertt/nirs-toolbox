function quantizefd(this,eventdata)
%QUANTIZEFD   

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

if isempty(this.reffracdelay)
    return;
end

% Quantize the coefficients
this.privfracdelay = quantizefd(this.filterquantizer,this.reffracdelay);


% [EOF]
