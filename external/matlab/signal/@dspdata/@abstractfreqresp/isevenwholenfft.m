function iseven = isevenwholenfft(this,Nfft,w)
%ISEVENWHOLENFFT   True if the length of the "whole" frequency response is even.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

if nargin < 2,
    [Nfft,nchans] = size(this.Data);
    w = this.Frequencies;
end

% To determine if the "whole" Nfft is EVEN check to see if the frequency
% vector includes the Nyquist (pi or Fs/2).
if this.NormalizedFrequency,
    fnyq = pi; 
else
    fnyq = this.getfs/2;
end

lastpt    = w(end);
freqrange = lastpt-w(1);
halfDeltaF  = freqrange/(Nfft-1)/2;
dist2nonnyq = lastpt - (fnyq-halfDeltaF);  % Distance to pt before Nyquist.
dist2nyq    = lastpt - fnyq;               % Distance to Nyquist.

iseven = false;
if  abs(dist2nyq) < abs(dist2nonnyq), % Assume EVEN "whole" NFFT
    iseven = true;
end

% [EOF]
