function [num, den, z, p] = besself(n, Wn, ftype, anaflag) %#ok
%BESSELF  Bessel analog filter design.
%   [B,A] = BESSELF(N,Wo) designs an N'th order lowpass analog Bessel
%   filter and returns the filter coefficients in length N+1 vectors B and
%   A.  The frequency Wo is the frequency up to which the group delay is
%   approximately constant.
%
%   The larger the filter order, the better the filter's group delay will
%   approximate a constant up to the frequency Wo.
%
%   Note that Bessel filters are lowpass only.
%
%   When used with three left-hand arguments, as in [Z,P,K] = BESSELF(...),
%   the zeros and poles are returned in length N column vectors Z and P,
%   and the gain in scalar K.
%
%   When used with four left-hand arguments, as in [A,B,C,D] =
%   BESSELF(...), state-space matrices are returned.
%
%   % Example 1:
%   %   Design a fifth-order analog lowpass Bessel filter with an
%   %   approximate constant group delay up to 10,000 rad/s and plot
%   %   the frequency response of the filter using freqs.
%
%   [b,a] = besself(5,10000);
%   freqs(b,a)                   % Plot frequency response
%
%   % Example 2:
%   %   Design a 5-th order Bessel analog filter and convert it to a
%   %   digital filter.
%
%   Fs =100;                            % Sampling frequency
%   [z,p,k] = besself(5,1000);          % Bessel analog filter design
%   [zd,pd,kd]=bilinear(z,p,k,Fs);      % Analog to digital mapping
%   [sos] = zp2sos(zd,pd,kd);           % Convert to SOS form
%   fvtool(sos)                         % Visualize the digital filter
%
%   See also BESSELAP, BUTTER, CHEBY1, CHEBY2, FREQZ, FILTER.

%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     [1] T. W. Parks and C. S. Burrus, Digital Filter Design,
%         John Wiley & Sons, 1987, chapter 7, section 7.3.3.

btype = 1;
if (nargin >= 3)
    % band-stop or high-pass
    btype = 3;
end
if length(Wn) == 2
    btype = btype + 1;
end
% Cast to enforce precision rules
n = signal.internal.sigcasttofloat(n,'double','besself','N',...
  'allownumeric');
Wn = signal.internal.sigcasttofloat(Wn,'double','besself','Wo',...
  'allownumeric');

% step 1: get analog, pre-warped frequencies
u = Wn;

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

% step 3: Get N-th order Bessel analog lowpass prototype
[z,p,k] = besselap(n);

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
        zinf = ltipack.getTolerance('infzero',true);
        [z,k] = ltipack.sszero(a,b,c,d,[],zinf);
        num = k * poly(z);
        num = [zeros(1,length(den)-length(num))  num];
    end
end

