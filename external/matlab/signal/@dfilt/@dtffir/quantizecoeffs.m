function quantizecoeffs(h,eventData)
% Quantize coefficients


%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

if isempty(h.refnum)
    return;
end

% Quantize the coefficients
h.privnum = quantizecoeffs(h.filterquantizer,h.refnum);

setmaxprod(h.filterquantizer, h);
setmaxsum(h.filterquantizer, h);

% [EOF]
