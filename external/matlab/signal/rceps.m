function [xhat, yhat] = rceps(x)
%RCEPS Real cepstrum.
%   RCEPS(X) returns the real cepstrum of the real sequence X.
%
%   [XHAT, YHAT] = RCEPS(X) returns both the real cepstrum XHAT and
%   YHAT which is a unique minimum-phase sequence that has the same real
%   cepstrum as X.
%
%   EXAMPLE: Show that YHAT is a unique minimum-phase sequence that has the
%            % same real cepstrum as X.
%            y = [4 1 5]; % Non-minimum phase sequence
%            [xhat,yhat] = rceps(y);
%            xhat2= rceps(yhat);
%            [xhat' xhat2']
%
%   See also CCEPS, ICCEPS, HILBERT, FFT.

%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     [1] A.V. Oppenheim and R.W. Schafer, Digital Signal
%         Processing, Prentice-Hall, 1975.
%     [2] Programs for Digital Signal Processing, IEEE Press,
%         John Wiley & Sons, 1979, algorithm 7.2.

narginchk(1,1);

% Check for valid input signal
chkinput(x);

n = length(x);
fftxabs = abs(fft(x));
if any(fftxabs == 0)
    error(message('signal:rceps:ZeroInFFT'));
end
xhat = real(ifft(log(fftxabs)));
if nargout > 1
    odd = fix(rem(n,2));
    wn = [1; 2*ones((n+odd)/2-1,1) ; ones(1-rem(n,2),1); zeros((n+odd)/2-1,1)];
    % class(yhat) must be the same as that of xhat to ensure precision rules
    yhat = zeros(size(x),class(xhat));
    yhat(:) = real(ifft(exp(fft(wn.*xhat(:)))));
end


%--------------------------------------------------------------------------
function chkinput(x)
% Check for valid input signal

if isempty(x) || issparse(x) || ~isreal(x),
    error(message('signal:rceps:invalidInput'));
end

