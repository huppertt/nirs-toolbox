function L = setblocklength(Hd, L)
%SETBLOCKLENGTH Overloaded set on the blocklength property.
  
%   Author: R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Get filter order
M = nstates(Hd);

% Set the fft coeffs as a column
fftcoeffs = fft(Hd.Numerator,L+M); % Be careful of Hd.Numerator = scalar
% since fft will always produce a row! Don't do: Hd.fftcoeffs = fft(Hd.Numerator(:),L+M);
Hd.fftcoeffs = fftcoeffs(:);
