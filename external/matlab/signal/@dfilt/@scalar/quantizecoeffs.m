function quantizecoeffs(h,eventData)
% Quantize coefficients


%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

if isempty(h.refgain)
    return;
end

% Quantize the coefficients
h.privgain = quantizecoeffs(h.filterquantizer,h.refgain);

% [EOF]
