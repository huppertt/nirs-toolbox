function quantizecoeffs(h,eventData)
% Quantize coefficients

%   Author(s): V. Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.
%   Date: 2004/06/06 16:54:31 $

if isempty(h.refBeta)
    return;
end

% Quantize the coefficients
[Allpass1q,Allpass2q,Betaq] = quantizecoeffs(h.filterquantizer,h.refAllpass1,h.refAllpass2,h.refBeta);
h.Allpass1q = Allpass1q;
h.Allpass2q = Allpass2q;
h.privBeta = Betaq;

% [EOF]
