function quantizecoeffs(Hd,eventData)
% Quantize coefficients

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

if isempty(Hd.refladder)
    return;
end

q = Hd.filterquantizer;
[latq, conjlatq, ladq] = quantizecoeffs(q,Hd.reflattice,conj(Hd.reflattice),Hd.refladder);
Hd.privlattice = latq;
Hd.privconjlattice = conjlatq;
Hd.privladder = ladq;


