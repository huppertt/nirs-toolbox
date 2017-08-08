function n = nsectionsfft(Hd)
%NSECTIONSFFT Returns the number of sections of the FFT.
  
%   Author: V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

nfft = nstates(Hd) + Hd.BlockLength;
n = log2(nfft);

if round(n)~=n,
    error(message('signal:dfilt:fftfir:nsectionsfft:InvalidParam'));
end
