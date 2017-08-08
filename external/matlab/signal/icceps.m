function x = icceps(xhat,nd)
%ICCEPS Inverse complex cepstrum.
%   ICCEPS(XHAT,ND) returns the inverse complex cepstrum of the real
%   sequence XHAT, removing ND samples of delay. If XHAT was obtained with
%   CCEPS(X), then the amount of delay that was added to X was the element
%   of round(unwrap(angle(fft(x)))/pi) corresponding to pi radians.
%
%   EXAMPLE: Use ICCEPS to compute the inverse complex cepstrum.
%   x = 1:10;
%   [xh,nd] = cceps(x);
%   % Use the delay parameter with icceps to invert the complex cepstrum
%   icceps(xh,nd)
%
%   See also CCEPS, RCEPS, HILBERT, and FFT.

%   Author(s): T. Krauss, 2/22/96
%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     [1] A.V. Oppenheim and R.W. Schafer, Digital Signal
%         Processing, Prentice-Hall, 1975.

narginchk(1,2);

% Check for valid input signal
chkinput(xhat);

if nargin<2
    nd = 0;
end
% Cast to enforce precision rules
nd = signal.internal.sigcasttofloat(nd,'double','icceps','ND',...
  'allownumeric');

logh = fft(xhat);
h = exp(real(logh)+1i*rcwrap(imag(logh),nd));
x = real(ifft(h));

%--------------------------------------------------------------------------
function x = rcwrap(y,nd)
%RCWRAP Phase wrap utility used by ICCEPS.
%   RCWRAP(X,nd) adds phase corresponding to integer lag.

n = length(y);
nh = fix((n+1)/2);
x = y(:).' + pi*nd*(0:(n-1))/nh;
if size(y,2)==1
    x = x.';
end

%--------------------------------------------------------------------------
function chkinput(x)
% Check for valid input signal

if isempty(x) || issparse(x) || ~isreal(x),
    error(message('signal:icceps:invalidInput', 'X'));
end
