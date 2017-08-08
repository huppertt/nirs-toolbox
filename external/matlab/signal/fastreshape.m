function [h, w] = fastreshape(h, fs, unitcircle, nfft)
%FASTRESHAPE   Reshapes according to the frequency range.
%   FASTRESHAPE(HOBJ, H, FS, UC, NFFT)
%
%   This is a helper function intended to be used by object methods.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

switch unitcircle
    case 1,  % 0-pi
        h = h(1:nfft/2, :);  % Return half since h contains 2*nfft points.
        w = freqz_freqvec(nfft/2, fs, 2);
        
    case 2,  % 0-2pi
        w = freqz_freqvec(nfft, fs, 1);
        
    case 3,  % -p-pi
        h = fftshift(h);
        w = freqz_freqvec(nfft, fs, 3);
end
w=w(:);

% [EOF]
