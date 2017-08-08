function [num, den, z, p] = cheby1(n, r, Wn, varargin)
%CHEBY1 Chebyshev Type I digital and analog filter design.
%   [B,A] = CHEBY1(N,R,Wp) designs an Nth order lowpass digital Chebyshev
%   filter with R decibels of peak-to-peak ripple in the passband. CHEBY1
%   returns the filter coefficients in length N+1 vectors B (numerator) and
%   A (denominator). The passband-edge frequency Wp must be 0.0 < Wp < 1.0,
%   with 1.0 corresponding to half the sample rate.  Use R=0.5 as a
%   starting point, if you are unsure about choosing R.
%
%   If Wp is a two-element vector, Wp = [W1 W2], CHEBY1 returns an order 2N
%   bandpass filter with passband  W1 < W < W2. [B,A] =
%   CHEBY1(N,R,Wp,'high') designs a highpass filter. [B,A] =
%   CHEBY1(N,R,Wp,'low') designs a lowpass filter. [B,A] =
%   CHEBY1(N,R,Wp,'stop') is a bandstop filter if Wp = [W1 W2].
%
%   When used with three left-hand arguments, as in [Z,P,K] = CHEBY1(...),
%   the zeros and poles are returned in length N column vectors Z and P,
%   and the gain in scalar K.
%
%   When used with four left-hand arguments, as in [A,B,C,D] = CHEBY1(...),
%   state-space matrices are returned.
%
%   CHEBY1(N,R,Wp,'s'), CHEBY1(N,R,Wp,'high','s') and
%   CHEBY1(N,R,Wp,'stop','s') design analog Chebyshev Type I filters. In
%   this case, Wp is in [rad/s] and it can be greater than 1.0.
%
%   % Example 1:
%   %   For data sampled at 1000 Hz, design a 9th-order lowpass Chebyshev
%   %   Type I filter with 5 dB of ripple in the passband, and a passband
%   %   edge frequency of 300Hz.
%
%   Wn = 300/500;               % Normalized passband edge frequency
%   [z,p,k] = cheby1(9,5,Wn);
%   [sos] = zp2sos(z,p,k);      % Convert to SOS form
%   h = fvtool(sos)             % Plot magnitude response
%
%   % Example 2:
%   %   Design a 2nd-order Chebyshev Type I band-pass filter which passes
%   %   frequencies between 0.2 and 0.5 with 3 dB of ripple in the
%   %   passband.
%
%   [b,a]=cheby1(2,3,[.2,.5]);      % Bandpass digital filter design
%   h = fvtool(b,a);                % Visualize filter
%
%   See also CHEB1ORD, CHEBY2, BUTTER, ELLIP, FREQZ, FILTER, DESIGNFILT.

%   Author(s): J.N. Little, 1-14-87
%   	   J.N. Little, 1-13-88, revised
%   	   L. Shure, 4-29-88, revised
%   	   T. Krauss, 3-25-93, revised

%   Copyright 1988-2013 The MathWorks, Inc.


%   References:
%     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%         John Wiley & Sons, 1987, chapter 7, section 7.3.3.

% Cast to enforce precision rules
n = signal.internal.sigcasttofloat(n,'double','cheby1','N',...
  'allownumeric');
r = signal.internal.sigcasttofloat(r,'double','cheby1','R',...
  'allownumeric');
Wn = signal.internal.sigcasttofloat(Wn,'double','Wp',...
  'allownumeric');

[btype,analog,errStr,msgobj] = iirchk(Wn,varargin{:});
if ~isempty(errStr), error(msgobj); end

if n>500
    error(message('signal:cheby1:InvalidRange'))
end

% step 1: get analog, pre-warped frequencies
if ~analog,
    fs = 2;
    u = 2*fs*tan(pi*Wn/fs);
else
    u = Wn;
end

Bw = [];
% step 2: convert to low-pass prototype estimate
if btype == 1	% lowpass
    Wn = u;
elseif btype == 2	% bandpass
    Bw = u(2) - u(1);
    Wn = sqrt(u(1)*u(2));	% center frequency
elseif btype == 3	% highpass
    Wn = u;
elseif btype == 4	% bandstop
    Bw = u(2) - u(1);
    Wn = sqrt(u(1)*u(2));	% center frequency
end

% step 3: Get N-th order Chebyshev type-I lowpass analog prototype
[z,p,k] = cheb1ap(n, r);

% Transform to state-space
[a,b,c,d] = zp2ss(z,p,k);

% step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
if btype == 1		% Lowpass
    [a,b,c,d] = lp2lp(a,b,c,d,Wn);
    
elseif btype == 2	% Bandpass
    [a,b,c,d] = lp2bp(a,b,c,d,Wn,Bw);
    
elseif btype == 3	% Highpass
    [a,b,c,d] = lp2hp(a,b,c,d,Wn);
    
elseif btype == 4	% Bandstop
    [a,b,c,d] = lp2bs(a,b,c,d,Wn,Bw);
end

% step5: Use Bilinear transformation to find discrete equivalent:
if ~analog,
    [a,b,c,d] = bilinear(a,b,c,d,fs);
end

if nargout == 4
    num = a;
    den = b;
    z = c;
    p = d;
else	% nargout <= 3
    % Transform to zero-pole-gain and polynomial forms:
    if nargout == 3
        [z,p,k] = ss2zp(a,b,c,d,1); %#ok
        z = cheb1zeros(btype,n,Wn,analog);
        num = z;
        den = p;
        z = k;
    else % nargout <= 2
        den = poly(a);
        num = cheb1num(btype,n,Wn,Bw,analog,den,r);
        % num = poly(a-b*c)+(d-1)*den;
    end
end

%---------------------------------
function b = cheb1num(btype,n,Wn,Bw,analog,den,Rp)
% This internal function returns more exact numerator vectors
% for the num/den case.
% normalize so |H(w)| == 0 or -Rp decibels:
if ~rem(n,2)
    g = 10^(-Rp/20);
else
    g = 1;
end
if analog
    switch btype
        case 1  % lowpass
            b = [zeros(1,n) n^(-n)];
            b = real( g*b*polyval(den,-1i*0)/polyval(b,-1i*0) );
        case 2  % bandpass
            b = [zeros(1,n) Bw^n zeros(1,n)];
            b = real( g*b*polyval(den,-1i*Wn)/polyval(b,-1i*Wn) );
        case 3  % highpass
            b = [1 zeros(1,n)];
            b = real( g*b*den(1)/b(1) );
        case 4  % bandstop
            r = 1i*Wn*((-1).^(0:2*n-1)');
            b = poly(r);
            b = real( g*b*polyval(den,-1i*0)/polyval(b,-1i*0) );
    end
else
    Wn = 2*atan2(Wn,4);
    switch btype
        case 1  % lowpass
            r = -ones(n,1);
            w = 0;
        case 2  % bandpass
            r = [ones(n,1); -ones(n,1)];
            w = Wn;
        case 3  % highpass
            r = ones(n,1);
            w = pi;
        case 4  % bandstop
            r = exp(1i*Wn*( (-1).^(0:2*n-1)' ));
            w = 0;
    end
    b = poly(r);
    kern = exp(-1i*w*(0:length(b)-1));
    b = real(g*b*(kern*den(:))/(kern*b(:)));
end

function z = cheb1zeros(btype,n,Wn,analog)
% This internal function returns more exact zeros.
if analog
    % for lowpass and bandpass, don't include zeros at +Inf or -Inf
    switch btype
        case 1  % lowpass
            z = zeros(0,1);
        case 2  % bandpass
            z = zeros(n,1);
        case 3  % highpass
            z = zeros(n,1);
        case 4  % bandstop
            z = 1i*Wn*((-1).^(0:2*n-1)');
    end
else
    Wn = 2*atan2(Wn,4);
    switch btype
        case 1  % lowpass
            z = -ones(n,1);
        case 2  % bandpass
            z = [ones(n,1); -ones(n,1)];
        case 3  % highpass
            z = ones(n,1);
        case 4  % bandstop
            z = exp(1i*Wn*( (-1).^(0:2*n-1)' ));
    end
end


