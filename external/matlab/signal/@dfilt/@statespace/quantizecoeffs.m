function quantizecoeffs(h,eventData)
% Quantize coefficients

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.
%   Date: 2004/04/12 23:59:35 $


% Quantize the coefficients
[Aq,Bq,Cq,Dq] = quantizecoeffs(h.filterquantizer,h.refA,h.refB,h.refC,h.refD);
h.privA = Aq;
h.privB = Bq;
h.privC = Cq;
h.privD = Dq;


