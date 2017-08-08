function [xhat,nd,xhat1] = cceps(x,n)
%CCEPS Complex cepstrum.
%   CCEPS(X) returns the complex cepstrum of the real sequence X using
%   Fourier transform as in its definition [1]. The input is altered to
%   have no phase discontinuity at +/- pi radians by the application of a
%   linear phase term (that is, it is circularly shifted (after
%   zero-padding) by some samples if necessary to have zero phase at pi
%   radians).  The program follows the procedure used in [2]
%
%   [XHAT,ND] = CCEPS(X) returns the number of samples ND of (circular)
%   delay added to X prior to finding the complex cepstrum.
%
%   [XHAT,ND,XHAT1] = CCEPS(X) returns the complex cepstrum in XHAT1 using 
%   an alternate factorization algorithm developed in [3]. The first
%   element in the data sequence cannot be zero.
%
%   CCEPS(X,N) zero-pads X to length N, and returns the length N complex 
%   cepstrum of X.
%
%   The comparison between the two algorithms can be seen from the
%   following table: 
%
%   ----------------------------------------------------------------------
%   |  Algorithm   |        Pros          |        Cons                  |
%   ----------------------------------------------------------------------
%   |  Fourier     |   Can be applied to  |  Needs phase unwrapping.     |
%   |              |   any signal.        |  Output is aliased.          |
%   ----------------------------------------------------------------------
%   |Factorization |  Does not need phase | The input signal has to have |
%   |              |  unwrapping.         | an all-zero Z transform.     |
%   |              |  No aliasing.        | Z transform should have no   |
%   |              |                      | zeros on the unit circle.    |
%   ----------------------------------------------------------------------
%
%   Generally, the results of the two algorithms cannot be used to verify
%   each other.  
%
%   The two results can be used to verify each other only when the first
%   element of data sequence is positive, the Z transform of the data
%   sequence has only zeros, all these zeros are inside the unit circle,
%   and the data sequence is long or padded with enough zeros.
%
%   EXAMPLE: Use CCEPS to show an echo.
%   Fs = 100; t = 0:1/Fs:1.27;
%   s1 = sin(2*pi*45*t); % 45Hz sine wave sampled at 100Hz
% 
%   % Add an echo of the signal, with half the amplitude, 0.2 seconds after 
%   % the beginning of the signal. 
%   s2 = s1 + 0.5*[zeros(1,20) s1(1:108)];
%
%   c = cceps(s2);
%   plot(t,c)
%   title('The peak at 0.2 shows the echo');
%
%   See also ICCEPS, RCEPS, HILBERT, and FFT.

%   Copyright 1988-2013 The MathWorks, Inc.

%   References: 
%     [1] Oppenheim, A.V. and Schafer, R.W.  Discrete-Time Signal 
%         Processing, Prentice-Hall, 1989.
%     [2] Programs for Digital Signal Processing, IEEE Press,
%         John Wiley & Sons, 1979, algorithm 7.1.
%     [3] Steiglitz, K. and Dickinson, B.  "Computation of the complex
%         cepstrum by factorization of the Z-transform," Proc. Int. Conf.
%         ASSP, 1977, 723-726

narginchk(1,2);

% Check for valid input signal
chkinput(x);

if nargin < 2
    h = fft(x);
else
    n = signal.internal.sigcasttofloat(n,'double','cceps','N',...
      'allownumeric');
    h = fft(x,n);
end
[ah,nd] = rcunwrap(angle(h));
logh = log(abs(h))+1i*ah;
xhat = real(ifft(logh));

if nargout==3
% use alternate rooting algorithm to check original result
% Oppenheim & Schafer, p.795 [1] and Steiglitz & Dickinson [3] 
%   - doesn't work with zeros on the unit circle
%   - only works for finite sequences that can be rooted 
    x = x(:).'; % Force vector to be a row
    r = roots(x);
    r = r(r~=0);
    
    if ~isempty(r(abs(r) == 1))
        error(message('signal:cceps:SignalErr'));
    end
    
    a = r(abs(r)<1);
    b = r(abs(r)>1);
    
    A = x(1)*prod(b);
    b = 1./b;
    
    if nargin == 2
      datalen = n;
    else
      datalen = length(x);
    end
    % class(xhat1) should be same as class(xhat) to enforce precision rules
    xhat1 = zeros(datalen,1,class(xhat)); %#ok<*ZEROLIKE>
    xhat1(1) = log(abs(A));
    
    % n > 0 contributions are from zeros inside the unit circle
    if ~isempty(a)
        n = 1:datalen-1;
        [a,n]=meshgrid(a,n');
        xhat1 = xhat1 + [0; -((a.^n)./n)*ones(size(a,2),1)];
    end

    % n < 0 contributions are from zeros outside the unit circle
    if ~isempty(b)
        n = -(datalen-1):-1;
        [b,n]=meshgrid(b,n');
        xhat1 = xhat1 + [0; ((b.^(-n))./n)*ones(size(b,2),1)];
    end

    if size(x,1)==1
        xhat1 = xhat1';
    end
end

%--------------------------------------------------------------------------
function [y,nd] = rcunwrap(x)
%RCUNWRAP Phase unwrap utility used by CCEPS.
%   RCUNWRAP(X) unwraps the phase and removes phase corresponding
%   to integer lag.  See also: UNWRAP, CCEPS.

n = length(x);
y = unwrap(x);
nh = fix((n+1)/2);

idx = nh+1; 
% Special case the index for scalar input.
if length(y) == 1
  idx = 1;
end
nd = round(y(idx)/pi);
y(:) = y(:)' - pi*nd*(0:(n-1))/nh;
% Cast to 'double' to enforce precision rules
nd = double(nd); % output nd has no bearing on single precision rules
%--------------------------------------------------------------------------
function chkinput(x)
% Check for valid input signal

signal.internal.sigcheckfloattype(x,'','cceps','X');
if isempty(x) || issparse(x) || ~isreal(x),
    error(message('signal:cceps:invalidInput', 'X'));
end

