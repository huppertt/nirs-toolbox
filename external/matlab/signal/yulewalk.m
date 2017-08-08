function [B,A] = yulewalk(na, ff, aa, npt, lap)
%YULEWALK Recursive filter design using a least-squares method.
%   [B,A] = YULEWALK(N,F,M) finds the N-th order recursive filter
%   coefficients B and A such that the filter:
%   	                      -1             -(n-1) 
%   	   B(z)   b(1) + b(2)z + .... + b(n)z
%   	   ---- = ---------------------------
%   	                      -1             -(n-1)
%   	   A(z)    1   + a(1)z + .... + a(n)z
%
%   matches the magnitude frequency response given by vectors F and M.
%   Vectors F and M specify the frequency and magnitude breakpoints for
%   the filter such that PLOT(F,M) would show a plot of the desired
%   frequency response. The frequencies in F must be between 0.0 and 1.0,
%   with 1.0 corresponding to half the sample rate. They must be in
%   increasing order and start with 0.0 and end with 1.0. 
%
%   % Example:
%   %   Design an 8th-order lowpass filter and overplot the desired  
%   %   frequency response with the actual frequency response.
%
%   f = [0 0.6 0.6 1];      % Frequency breakpoints 
%   m = [1 1 0 0];          % Magnitude breakpoints
%   [b,a] = yulewalk(8,f,m);% Filter design using a least-squares method
%   [h,w] = freqz(b,a,128); % Frequency response of filter
%   plot(f,m,w/pi,abs(h),'--')
%   legend('Ideal','yulewalk Designed')
%   title('Comparison of Frequency Response Magnitudes')
%
%   See also FIR1, BUTTER, CHEBY1, CHEBY2, ELLIP, FREQZ and FILTER.

%   The YULEWALK function performs a least squares fit in the time
%   domain. The denominator coefficients {a(1),...,a(NA)} are computed
%   by the so called "modified Yule Walker" equations, using NR
%   correlation coefficients computed by inverse Fourier transformation
%   of the specified frequency response H.
%   The numerator is computed by a four step procedure. First, a numerator
%   polynomial corresponding to an additive decomposition of the power 
%   frequency response is computed. Next, the complete frequency response
%   corresponding to the numerator and denominator polynomials is
%   evaluated. Then a spectral factorization technique is used to
%   obtain the impulse response of the filter. Finally, the numerator
%   polynomial is obtained by a least squares fit to this impulse
%   response. For a more detailed explanation of the algorithm see 
%   B. Friedlander and B. Porat, "The Modified Yule-Walker Method
%   of ARMA Spectral Estimation," IEEE Transactions on Aerospace
%   Electronic Systems, Vol. AES-20, No. 2, pp. 158-173, March 1984.
%
%   Perhaps the best way to get familiar with the proper use of YULEWALK
%   is to examine carefully the demonstration file FILTDEMO.M, which
%   provides a detailed example.
%
%   See also BUTTER, CHEBY1, CHEBY2, ELLIP, FIR2, FIRLS, MAXFLAT and
%   FIRPM.

%   Author(s): B. Friedlander, 7-16-85
%   	   L. Shure, 3-5-87, modified
%   	   C. Denham, 7/26/90, repaired
%   Copyright 1988-2013 The MathWorks, Inc.
%     

if (nargin < 3 || nargin > 5)
   error(message('signal:yulewalk:Nargchk'))
end
% Cast to enforce Precision Rules
na = signal.internal.sigcasttofloat(na,'double','yulewalk','N',...
  'allownumeric');

if (nargin > 3)
   npt = signal.internal.sigcasttofloat(npt,'double','yulewalk','NPT',...
     'allownumeric');
   if round(2 .^ round(log(npt)/log(2))) ~= npt
	% NPT is not an even power of two
      npt = round(2.^ceil(log(npt)/log(2)));
   end 
end
if (nargin < 4)
   npt = 512;
end
if (nargin < 5)
   lap = fix(npt/25);
end
% Checks if LAP is a valid numeric input
lap = signal.internal.sigcasttofloat(lap,'double','yulewalk',...
  'LAP','allownumeric');

% Checks if M  and F are valid
signal.internal.sigcheckfloattype(ff,'','yulewalk','F');
signal.internal.sigcheckfloattype(aa,'','yulewalk','M');

[mf,nf] = size(ff);
[mm,nn] = size(aa);
if mm ~= mf || nn ~= nf
   error(message('signal:yulewalk:InvalidDimensions'))
end
nbrk = max(mf,nf);
if mf < nf
   ff = ff';
   aa = aa';
end

if abs(ff(1)) > eps || abs(ff(nbrk) - 1) > eps
   error(message('signal:yulewalk:InvalidFreqVec'))
end

% interpolate breakpoints onto large grid

npt = npt + 1;  % For [dc 1 2 ... nyquist].
Ht = zeros(1,npt);

nint=nbrk-1;
df = diff(ff'); 
if (any(df < 0))
   error(message('signal:yulewalk:MonotonicFreq'))
end

nb = 1;
Ht(1)=aa(1);
for i=1:nint
    if df(i) == 0
       nb = nb - lap/2;
       ne = nb + lap;
    else
       ne = fix(ff(i+1)*npt);
    end
    if (nb < 0 || ne > npt)
       error(message('signal:yulewalk:SignalErr'))
    end
    j=nb:ne;
    if ne == nb
        inc = 0;
    else
        inc = (j-nb)/(ne-nb);
    end
    Ht(nb:ne) = inc*aa(i+1) + (1 - inc)*aa(i);
    nb = ne + 1;
end
Ht = [Ht Ht(npt-1:-1:2)];
n = length(Ht);
n2 = fix((n+1)/2);
nb = na;
nr = 4*na;
nt = 0:1:nr-1;

% compute correlation function of magnitude squared response

R = real(ifft(Ht .* Ht));

R  = R(1:nr).*(0.54+0.46*cos(pi*nt/(nr-1)));     % pick NR correlations 

% Form window to be used in extracting the right "wing" of two-sided
% covariance sequence.
Rwindow = [1/2 ones(1,n2-1) zeros(1,n-n2)]; 

A = polystab(denf(R,na));            	% compute denominator

Qh = numf([R(1)/2,R(2:nr)],A,na);	% compute additive decomposition

Ss = 2*real(freqz(Qh,A,n,'whole'))';    % compute impulse response
hh = ifft(exp(fft(Rwindow.*ifft(log(Ss)))));
B  = real(numf(hh(1:nr),A,nb));

function b = numf(h,a,nb)
%NUMF	Find numerator B given impulse-response h of B/A and denominator A
%   NB is the numerator order.  This function is used by YULEWALK.
  
%       Was Revision: 1.3, Date: 1994/01/25 17:59:33 

nh = max(size(h)); 
impr = filter(1,a,[1 zeros(1,nh-1)]);
b = h/toeplitz(impr,[1 zeros(1,nb)])';

function A = denf(R,na)
%DENF	Compute denominator from covariances.
%   A = DENF(R,NA) computes order NA denominator A from covariances 
%   R(0)...R(nr) using the Modified Yule-Walker method.  
%   This function is used by YULEWALK.

%       Was Revision: 1.5, Date: 1994/01/25 17:59:00

nr = max(size(R));
Rm = toeplitz(R(na+1:nr-1),R(na+1:-1:2));
Rhs = - R(na+2:nr);
A = [1 Rhs/Rm'];

