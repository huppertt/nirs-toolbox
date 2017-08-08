function n = filtord(b,varargin)
%FILTORD Filter order
%   N = FILTORD(B,A) returns the order, N, of the filter:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   given numerator and denominator coefficients in vectors B and A. 
%
%   N = FILTORD(SOS) returns order, N, of the filter specified using the
%   second order sections matrix SOS. SOS is a Kx6 matrix, where the number
%   of sections, K, must be greater than or equal to 2. Each row of SOS
%   corresponds to the coefficients of a second order filter. From the
%   transfer function displayed above, the ith row of the SOS matrix
%   corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   N = FILTORD(D) returns order, N, of the digital filter D. You design a
%   digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   % Example 1:
%   %   Create a 25th-order lowpass FIR filter and verify its order using 
%   %   filtord.
%
%   b = fir1(25,0.45);
%   n = filtord(b)    
%
%   % Example 2:
%   %   Create a 6th-order lowpass IIR filter using second order sections
%   %   and verify its order using filtord.
%
%   [z,p,k] = butter(6,0.35);
%   SOS = zp2sos(z,p,k);	% second order sections matrix
%   n = filtord(SOS)	    
%
%   % Example 3:
%   %   Use the designfilt function to design a highpass IIR digital filter 
%   %   with order 8, passband frequency of 75 KHz, and a passband ripple 
%   %   of 0.2 dB. Use filtord to get the filter order. 
%  
%   D = designfilt('highpassiir', 'FilterOrder', 8, ...
%            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%            'SampleRate', 200e3);
%
%   n = filtord(D)
%
%   See also FVTOOL, ISALLPASS, ISLINPHASE, ISMAXPHASE, ISMINPHASE, ISSTABLE

%   Copyright 2012-2013 The MathWorks, Inc.

narginchk(1,2);

a = 1; % Assume FIR for now

if all(size(b)>[1 1])
  % Input is a matrix, check if it is a valid SOS matrix
  if size(b,2) ~= 6
    error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
  end
  if nargin > 1
    error(message('signal:signalanalysisbase:invalidNumInputs'));
  end
  
  [b,a] = sos2tf(b); % get transfer function
  
elseif ~isempty(varargin)
  
  a = varargin{1};
  
  if all(size(a)>[1 1])
    error(message('signal:signalanalysisbase:inputnotsupported'));
  end
end

if ~isempty(b) && max(abs(b))~=0
  b = b/max(abs(b));
end
if ~isempty(a) && max(abs(a))~=0
  a = a/max(abs(a));
end

% Remove trailing "zeros" of a & b.
if ~isempty(b)
  b = b(1:find(b~=0, 1, 'last' ));
end
if ~isempty(a)
  a = a(1:find(a~=0, 1, 'last' ));
end

n = max(length(b),length(a)) - 1;
