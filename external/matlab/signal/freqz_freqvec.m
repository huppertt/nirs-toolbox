function w = freqz_freqvec(nfft, Fs, s)
%FREQZ_FREQVEC Frequency vector for calculating filter responses.
%   This is a helper function intended to be used by FREQZ.
%
%   Inputs:
%       nfft    -   The number of points
%       Fs      -   The sampling frequency of the filter
%       s       -   1 = 0-2pi, 2 = 0-pi, 3 = -pi-pi

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 2,  Fs = []; end
if nargin < 3,  s  = 2; end
if isempty(Fs), Fs = 2*pi; end

switch s
    case 1,  % 0-2pi
        deltaF = Fs/nfft;
        w = linspace(0,Fs-deltaF,nfft);
        
        % There can still be some minor round off errors.  Fix the known points,
        % those near pi and 2pi.
        if rem(nfft, 2)
            w((nfft+1)/2) = Fs/2-Fs/(2*nfft);
            w((nfft+1)/2+1) = Fs/2+Fs/(2*nfft);
        else
            % Make sure we hit Fs/2 exactly for the 1/2 nyquist point.
            w(nfft/2+1) = Fs/2;
        end
        w(nfft) = Fs-Fs/nfft;

    case 2,  % 0-pi
        deltaF = Fs/2/nfft;
        w = linspace(0,Fs/2-deltaF,nfft);
        
        w(nfft) = Fs/2-Fs/2/nfft;
        
    case 3, % -pi-pi
        deltaF = Fs/nfft;
        
        if rem(nfft,2), % ODD, don't include Nyquist.
            wmin = -(Fs - deltaF)/2;
            wmax = (Fs - deltaF)/2;
            
        else            % EVEN include Nyquist point in the negative freq.
            wmin = -Fs/2;
            wmax = Fs/2 - deltaF;
        end
        w = linspace(wmin, wmax, nfft);
        if rem(nfft, 2) % ODD
            w((nfft+1)/2) = 0;
        else
            w(nfft/2+1) = 0;
        end
end

% [EOF]
