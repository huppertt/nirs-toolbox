function beta = kaiserBeta(atten)
%This function is for internal use only. It may be removed in the future.

%kaiserBeta     Beta parameter for Kaiser window
%   BETA = kaiserBeta(ATTEN) returns the beta parameter for a Kaiser window
%   based on the sidelobe attenuation specified in ATTEN (in dB).
%
%   % Example:
%   %   Create a length-10 Kaiser window which can achieve a sidelobe
%   %   attenuation of 30 dB.
%
%   beta = signal.internal.kaiserBeta(30);
%   w = kaiser(10,beta)

%   Copyright 2012 The MathWorks, Inc.

beta = 0.1102*(atten-8.7).*(atten>50) + ...
    (0.5842*(atten-21).^0.4 + ...
    0.07886*(atten-21)).*(atten>=21 & atten<=50);

% [EOF]
