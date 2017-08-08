function quantizecoeffs(this,eventData)
%QUANTIZECOEFFS   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

if isempty(this.refallpasscoeffs),
    return;
end

% Quantize the coefficients
this.privallpasscoeffs = quantizecoeffs(this.filterquantizer,this.refallpasscoeffs);

% [EOF]
