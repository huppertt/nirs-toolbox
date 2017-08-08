function N = impzlength(b,varargin)
%IMPZLENGTH Impulse response length of digital filter
%   L = IMPZLENGTH(B,A) returns the length, L, of the impulse response of
%   the filter:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   given numerator and denominator coefficients in vectors B and A.
%
%   L = IMPZLENGTH(SOS) returns the approximate length, L, of the impulse
%   response of the filter specified using a second order sections matrix
%   SOS. SOS is a Kx6 matrix, where the number of sections, K, must be
%   greater than or equal to 2. Each row of SOS corresponds to the
%   coefficients of a second order filter. From the transfer function
%   displayed above, the ith row of the SOS matrix corresponds to [bi(1)
%   bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   L = IMPZLENGTH(D) returns the approximate length, L, of the impulse
%   response of the digital filter, D. You design a digital filter, D, by
%   calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   L = IMPZLENGTH(...,TOL) specifies a tolerance for greater or less
%   accuracy.  By default, TOL = 5e-5.
%
%   % Example 1:
%   %   Design a lowpass FIR filter with normalized cut-off frequency at
%   %   0.3 and determine the length of its impulse response.
%
%   b=fircls1(54,0.3,0.02,0.008);
%   impzlength(b)                       % length of impulse response
%
%   % Example 2:
%   %   Design a 5th order lowpass elliptic IIR filter and determine the
%   %   length of its impulse response.
%
%   [b,a] = ellip(5,0.5,20,0.4);
%   impzlength(b,a)                     % length of impulse response
%
%   % Example 3:
%   %   Design a Butterworth highpass IIR filter, represent its coefficients
%   %   using second order sections, and determine the length of its
%   %   impulse response.
%
%   [z,p,k] = butter(6,0.7,'high');
%   SOS = zp2sos(z,p,k);
%   impzlength(SOS)                     % length of impulse response
%
%   % Example 4:
%   %   Use the designfilt function to design a highpass IIR digital filter 
%   %   with order 8, passband frequency of 75 KHz, and a passband ripple 
%   %   of 0.2 dB. Sample rate is 200 KHz. Determine the length of its
%   %   impulse response.
%  
%   D = designfilt('highpassiir', 'FilterOrder', 8, ...
%            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%            'SampleRate', 200e3);
%
%   impzlength(D)
%
%   See also IMPZ.

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,3)

isTF = true; % True if dealing with a transfer function

if all(size(b)>[1 1])
    % Input is a matrix, check if it is a valid SOS matrix
    if size(b,2) ~= 6
        error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
    end
    isTF = false; % SOS instead of transfer function
    
    if nargin > 1
        tol = varargin{1};
    else
        tol = .00005;
    end
    % Checks if SOS is a valid numeric data input
    signal.internal.sigcheckfloattype(b,'','impzlength','SOS');
    % Cast to enforce precision rules
    tol = signal.internal.sigcasttofloat(tol,'double','impzlength',...
      'tol','allownumeric');
end

if isTF
    if nargin == 1
        a = 1;
    else
        a = varargin{1};
    end
    
    % Checks if B and A are valid numeric data inputs
    signal.internal.sigcheckfloattype(b,'','impzlength','B');
    signal.internal.sigcheckfloattype(a,'','impzlength','A');

    if nargin < 3
        tol = .00005;
    else
        tol = varargin{2};
    end
    % Cast to enforce precision rules
    tol = signal.internal.sigcasttofloat(tol,'double','impzlength',...
      'tol','allownumeric');
    
    % Determine if filter is FIR
    if signalpolyutils('isfir',b,a),
        N = length(b);
    else
        indx= find(b, 1);
        if isempty(indx),
            delay = 0;
        else
            delay=indx-1;
        end
        p = roots(a);
        if any(abs(p)>1.0001),
            N = unstable_length(p);
        else
            N = stableNmarginal_length(p,tol,delay);
        end
        % Cast to enforce precision rules
        N = double(N);
        N = max(length(a)+length(b)-1,N);
        
        % Always return an integer length
        N = floor(N);
    end
else
    N = lclsosimpzlength(b,tol);
    % Cast to enforce precision rules
    N = double(N);
end
%-------------------------------------------------------------------------
function N = unstable_length(p)
% Determine the length for an unstable filter
ind = abs(p)>1;
N = 6/log10(max(abs(p(ind))));% 1000000 times original amplitude


%-------------------------------------------------------------------------
function N = stableNmarginal_length(p,tol,delay)
% Determine the length for an unstable filter

%minimum height is .00005 original amplitude:
ind = find(abs(p-1)<1e-5);
p(ind) = -p(ind);    % treat constant as Nyquist
ind = find(abs(abs(p)-1)<1e-5);
periods = 5*max(2*pi./abs(angle(p(ind)))); % five periods
p(ind) = [];   % get rid of unit circle poles
[maxp,maxind] = max(abs(p));
if isempty(p)   % pure oscillator
    N = periods;
elseif isempty(ind)   % no oscillation
    N = mltplcty(p,maxind)*log10(tol)/log10(maxp) + delay;
else    % some of both
    N = max(periods, ...
        mltplcty(p,maxind)*log10(tol)/log10(maxp) ) + delay;
end

%-------------------------------------------------------------------------
function m = mltplcty( p, ind, tol)
%MLTPLCTY  Multiplicity of a pole
%   MLTPLCTY(P,IND,TOL) finds the multiplicity of P(IND) in the vector P
%   with a tolerance of TOL.  TOL defaults to .001.

if nargin<3
    tol = .001;
end

[mults,indx]=mpoles(p,tol);

m = mults(indx(ind));
for i=indx(ind)+1:length(mults)
    if mults(i)>m
        m = m + 1;
    else
        break;
    end
end

%--------------------------------------------------------------------------
function len = lclsosimpzlength(sos,tol)

% Initialize length
firlen=1;
iirlen=1;

% Convert the filter to a transfer function.
for k=1:size(sos,1)
    
    % Get the transfer function coefficients
    b=sos(k,1:3);
    a=sos(k,4:6);
    
    if signalpolyutils('isfir',b,a),
        % Add the length of each FIR section
        firlen = firlen + length(b) - 1;
    else
        % Keep the maximum length of all IIR sections
        iirlen = max(iirlen, impzlength(b,a,tol));
    end
end

% Use the longest of FIR or IIR
len=max(firlen,iirlen);

