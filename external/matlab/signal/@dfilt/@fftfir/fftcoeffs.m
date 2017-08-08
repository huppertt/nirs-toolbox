function c = fftcoeffs(Hd)
%FFTCOEFFS  Get the FFT coefficients used for the filtering.
%   C = FFTCOEFFS(Hd) returns the frequency-domain coefficients used in
%   the filtering. These coefficients are computed and stored prior to
%   filtering.
  
%   Author: R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Always return a row to be consistent with h.Numerator
% Use the first columns as we may have multipel columns repeated for
% filtering purposes
c = Hd.fftcoeffs(:,1).';
