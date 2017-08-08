function [H,W] = ispectrumshift(this,H,W)
%ISPECTRUMSHIFT   Inverse of SPECTRUMSHIFT.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

if nargin == 1,
    H = this.Data;
    W = this.Frequencies;
end

[nfft,nchans] = size(H);

% Determine half the number of FFT points.
if rem(nfft,2),
    halfNfft = (nfft+1)/2;  % ODD
    negEndPt = halfNfft;
    halfDeltaF = max(abs(diff(W)))/2;  % There's half of delta F on both sides of Nyquist.
else
    halfNfft = (nfft/2)+1;  % EVEN
    negEndPt = halfNfft-1;
    halfDeltaF = 0;  % Nyquist point is included.

    % Move the Nyquist point to the left-hand side (neg freq) as expected
    % by ifftshift.
    H = [H(end,:); H(1:end-1,:)];
end

% Convert to plot + frequencies only.
H = ifftshift(H);  

if this.normalizedFrequency,   Fn = pi;
else                           Fn = getfs(this)/2;
end
W = [W(negEndPt:end); -flipud(W(1:negEndPt-1))+Fn-halfDeltaF];

if nargout == 0,
    this.Data = H;
    this.Frequencies = W;
end

% [EOF]
