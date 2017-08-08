function [num, den, z, p] = cheby2(n, r, Wn, varargin)
%CHEBY2 Chebyshev Type II digital and analog filter design.
%   [B,A] = CHEBY2(N,R,Wst) designs an Nth order lowpass digital Chebyshev
%   filter with the stopband ripple R decibels down and stopband-edge
%   frequency Wst.  CHEBY2 returns the filter coefficients in length N+1
%   vectors B (numerator) and A (denominator). The stopband-edge frequency
%   Wst must be 0.0 < Wst < 1.0, with 1.0 corresponding to half the sample
%   rate.  Use R = 20 as a starting point, if you are unsure about choosing
%   R.
%
%   If Wst is a two-element vector, Wst = [W1 W2], CHEBY2 returns an order
%   2N bandpass filter with passband  W1 < W < W2. [B,A] =
%   CHEBY2(N,R,Wst,'high') designs a highpass filter. [B,A] =
%   CHEBY2(N,R,Wst,'low') designs a lowpass filter. [B,A] =
%   CHEBY2(N,R,Wst,'stop') is a bandstop filter if Wst = [W1 W2].
%
%   When used with three left-hand arguments, as in [Z,P,K] = CHEBY2(...),
%   the zeros and poles are returned in length N column vectors Z and P,
%   and the gain in scalar K.
%
%   When used with four left-hand arguments, as in [A,B,C,D] = CHEBY2(...),
%   state-space matrices are returned.
%
%   CHEBY2(N,R,Wst,'s'), CHEBY2(N,R,Wst,'high','s') and
%   CHEBY2(N,R,Wst,'stop','s') design analog Chebyshev Type II filters. In
%   this case, Wst is in [rad/s] and it can be greater than 1.0.
%
%   % Example 1:
%   %   For data sampled at 1000 Hz, design a ninth-order lowpass
%   %   Chebyshev Type II filter with stopband attenuation 40 dB down from
%   %   the passband and a stopband edge frequency of 300Hz.
%
%   Wn = 300/500;               % Normalized stopband edge frequency
%   [z,p,k] = cheby2(9,40,Wn);
%   [sos] = zp2sos(z,p,k);      % Convert to SOS form
%   h = fvtool(sos)             % Plot magnitude response
%
%   % Example 2:
%   %   Design a 6th-order Chebyshev Type II band-pass filter which passes
%   %   frequencies between 0.2 and 0.5 and with stopband attenuation 80 dB
%   %   down from the passband.
%
%   [b,a]=cheby2(6,80,[.2,.5]);     % Bandpass digital filter design
%   h = fvtool(b,a);                % Visualize filter
%
%   See also CHEB2ORD, CHEBY1, BUTTER, ELLIP, FREQZ, FILTER, DESIGNFILT.

%   Author(s): J.N. Little, 1-14-87
%   	   J.N. Little, 1-13-88, revised
%   	   L. Shure, 4-29-88, revised
%   	   T. Krauss, 3-25-93, revised
%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%         John Wiley & Sons, 1987, chapter 7, section 7.3.3.

% Cast to enforce precision rules
n = signal.internal.sigcasttofloat(n,'double','cheby2','N',...
  'allownumeric');
r = signal.internal.sigcasttofloat(r,'double','cheby2','R',...
  'allownumeric');
Wn = signal.internal.sigcasttofloat(Wn,'double','cheby2','Wst',...
  'allownumeric');

[btype,analog,errStr,msgobj] = iirchk(Wn,varargin{:});
if ~isempty(errStr), error(msgobj); end

if n>500
    error(message('signal:cheby2:InvalidRange'))
end

% step 1: get analog, pre-warped frequencies
if ~analog,
    fs = 2;
    u = 2*fs*tan(pi*Wn/fs);
else
    u = Wn;
end

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

% step 3: Get N-th order Chebyshev type-II lowpass analog prototype
[z,p,k] = cheb2ap(n, r);

% Transform to state-space
[a,b,c,d] = zp2ss(z,p,k);

% step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
if btype == 1		% Lowpass
    [a,b,c,d] = lp2lp(a,b,c,d,Wn);
    
elseif btype == 2	% Bandpass
    [a,b,c,d] = lp2bp(a,b,c,d,Wn,Bw);
    
elseif btype == 3	% Highpass
    [a,b,c,d] = lp2hp(a,b,c,d, Wn);
    
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
        [z,p,k] = ss2zp(a,b,c,d,1);
        num = z;
        den = p;
        z = k;
    else % nargout <= 2
        den = poly(a);
        if analog
            zeroLimitFlag = true;
        else
            zeroLimitFlag = false;
        end
        zinf = ltipack.getTolerance('infzero',zeroLimitFlag);
        [z,k] = ltipack.sszero(a,b,c,d,[],zinf);
        num = k * poly(z);
        num = [zeros(1,length(den)-length(num)) num];
    end
end

