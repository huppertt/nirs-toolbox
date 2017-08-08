function varargout = impz(b,varargin)
%IMPZ Impulse response of digital filter
%   [H,T] = IMPZ(B,A) computes the impulse response of the filter:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   given numerator and denominator coefficients in vectors B and A.
%
%   [H,T] = IMPZ(SOS) computes the impulse response of the filter specified
%   using the second order sections matrix SOS. SOS is a Kx6 matrix, where
%   the number of sections, K, must be greater than or equal to 2. Each row
%   of SOS corresponds to the coefficients of a second order filter. From
%   the transfer function displayed above, the ith row of the SOS matrix
%   corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   [H,T] = IMPZ(D) computes the impulse response of the digital filter, D.
%   You design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   In all cases the number of samples is chosen internally, and IMPZ
%   returns the response in column vector H and a vector of times (or
%   sample intervals) in T (T = [0 1 2 ...]').
%
%   [H,T] = IMPZ(...,N) computes N samples of the impulse response. If N is
%   a vector of integers, the impulse response is computed only at those
%   integer values (0 is the origin).
%
%   [H,T] = IMPZ(...,N,Fs) computes N samples and scales T so that samples
%   are spaced 1/Fs units apart.  Fs is 1 by default.
%
%   [H,T] = IMPZ(...,[],Fs) chooses the number of samples for you and
%   scales T so that samples are spaced 1/Fs units apart.
%
%   IMPZ(...) with no output arguments plots the impulse response using
%   STEM(T,H) in the current figure window.
%
%   % Example 1:
%   %   Design a lowpass FIR filter with normalized cut-off frequency at
%   %   0.3 and determine its impulse response.
%
%   b=fircls1(54,0.3,0.02,0.008);
%   impz(b)
%
%   % Example 2:
%   %   Design a 5th order lowpass elliptic IIR filter and determine its
%   %   impulse response.
%
%   [b,a] = ellip(5,0.5,20,0.4);
%   impz(b,a)
%
%   % Example 3:
%   %   Design a Butterworth highpass IIR filter, represent its coefficients
%   %   using second order sections, and display its impulse response.
%
%   [z,p,k] = butter(6,0.7,'high');
%   SOS = zp2sos(z,p,k);
%   impz(SOS)
%
%   % Example 4:
%   %   Use the designfilt function to design a highpass IIR digital filter 
%   %   with order 8, passband frequency of 75 KHz, and a passband ripple 
%   %   of 0.2 dB. Sample rate is 200 KHz. Visualize the impulse response 
%   %   using 256 samples.
%  
%   D = designfilt('highpassiir', 'FilterOrder', 8, ...
%            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%            'SampleRate', 200e3);
%
%   impz(D,256)
%
%   See also IMPULSE.
%
%   Note:  IMPULSE is in the Controls Toolbox for continuous systems.

%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(1,4)

isTF = true; % True if dealing with a transfer function

if all(size(b)>[1 1])
  % Input is a matrix, check if it is a valid SOS matrix
  if size(b,2) ~= 6
    error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
  end
  isTF = false; % SOS instead of transfer function
  
  if nargin > 1
    N = varargin{1};
  else
    N = [];
  end
  if nargin > 2
    Fs = varargin{2};
  else
    Fs = 1;
  end
  % Checks if SOS is a valid numeric data input
  signal.internal.sigcheckfloattype(b,'','impz','SOS');
  % Cast to enforce precision rules
  N = signal.internal.sigcasttofloat(N,'double','impz','N','allownumeric');
  Fs = signal.internal.sigcasttofloat(Fs,'double','impz','Fs',...
    'allownumeric');

end

% If input is a transfer function b,a -----
if isTF
  if nargin > 1
    a = varargin{1};
    if all(size(a)>[1 1])
      error(message('signal:signalanalysisbase:inputnotsupported'));
    end
  else
    a = 1;
  end
  
  % Checks if B and A are valid numeric data inputs
  signal.internal.sigcheckfloattype(b,'','impz','B');
  signal.internal.sigcheckfloattype(a,'','impz','A');

  if nargin > 2
    N = varargin{2};
  else
    N = [];
  end
  
  if nargin > 3
    Fs = varargin{3};
  else
    Fs = 1;
  end
  % Cast to enforce precision rules
  N = signal.internal.sigcasttofloat(N,'double','impz','N','allownumeric');
  Fs = signal.internal.sigcasttofloat(Fs,'double','impz','Fs',...
    'allownumeric');
end

% Compute time vector
M = 0;  NN = [];
if isempty(N)
  % if not specified, determine the length
  if isTF
    N = impzlength(b,a,.00005);
  else
    N  = impzlength(b,.00005);
  end
elseif length(N)>1  % vector of indices
  NN = round(N);
  N = max(NN)+1;
  M = min(min(NN),0);
end

t = (M:(N-1))';

% Filter an impulse with the input coefficients
x = double(t==0);
if isTF
  h = filter(b,a,x);
else
  h = sosfilt(b,x);
end

if ~isempty(NN),
  h = h(NN-M+1);
  t = t(NN-M+1);
end

t = t/Fs;

if nargout,
  % Cast to enforce precision rules
  if isa(h,'single')
    t = single(t);
  end
  varargout = {h, t};
else
  timezplot(t,h,Fs,'Impulse');
end
