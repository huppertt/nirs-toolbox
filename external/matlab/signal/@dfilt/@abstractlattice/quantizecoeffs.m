function quantizecoeffs(Hd,eventData)
%QUANTIZECOEFFS Quantize coefficients

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

% Quantize the coefficients
Hd.privlattice = quantizecoeffs(Hd.filterquantizer,Hd.reflattice);
Hd.privconjlattice = quantizecoeffs(Hd.filterquantizer,conj(Hd.reflattice));

