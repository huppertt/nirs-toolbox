function bw = enbw(window, fs)
%ENBW   Equivalent noise bandwidth
%   BW = ENBW(WINDOW) returns the two-sided equivalent noise bandwidth for
%   a uniformly sampled window whose coefficients are specified in the
%   vector WINDOW. This bandwidth is normalized to the noise power per
%   frequency bin.
%
%   BW = ENBW(WINDOW, Fs) returns the two-sided equivalent noise bandwidth
%   (in Hz) for a uniformly sampled window whose coefficients are specified
%   in the vector WINDOW, where Fs is the sampling rate of the window.
%
%   % Example 1:
%   %   Compute the equivalent noise bandwidth of a Hann window
%   bw1 = enbw(hann(10000))
%
%   % Example 2:
%   %   Compute the equivalent noise bandwidth (in Hz) of a Hann window
%   %   sampled at 44.1 kHz.
%   bw2 = enbw(hann(10000), 44.1e3)

%   Reference:
%     [1] fredric j. harris [sic], On the Use of Windows for Harmonic
%         Analysis with the Discrete Fourier Transform, Proceedings of
%         the IEEE, Vol. 66, No. 1, January 1978.  Eqn 11, 15.

%   Copyright 2012-2013 The MathWorks, Inc.

% two-sided ENBW computation of a window.
validateattributes(window,{'numeric'},{'real','vector'}, ...
    'enbw','WINDOW',1);

% compute normalized ENBW
bw = (rms(window)/mean(window))^2;

% if Fs is specified, scale by the frequency bin width
if nargin > 1
    validateattributes(fs,{'numeric'},{'real','positive','scalar'}, ...
        'enbw','Fs',2);
    bw = bw * double(fs) / length(window);
end

