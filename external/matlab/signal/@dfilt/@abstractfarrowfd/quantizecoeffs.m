function quantizecoeffs(h,eventData)
% Quantize coefficients


%   Author(s): R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.

if isempty(h.refcoeffs)
    return;
end

% Quantize the coefficients
h.privcoeffs = quantizecoeffs(h.privfilterquantizer,h.refcoeffs);

setmaxprod(h.filterquantizer, h);
setmaxsum(h.filterquantizer, h);

% [EOF]
