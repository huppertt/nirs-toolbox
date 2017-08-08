function f = islinphase(b,varargin)
%ISLINPHASE  True for linear phase filter
%   FLAG = ISLINPHASE(B,A) returns a logical output, FLAG, equal to TRUE if
%   the filter specified by numerator coefficients B, and denominator
%   coefficients A, is linear phase. Input vectors B, and A define a filter
%   with transfer function:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   FLAG = ISLINPHASE(SOS) returns TRUE if the filter specified using the
%   second order sections matrix, SOS, is linear phase. SOS is a Kx6
%   matrix, where the number of sections, K, must be greater than or equal
%   to 2. Each row of SOS corresponds to the coefficients of a second order
%   filter. From the transfer function displayed above, the ith row of the
%   SOS matrix corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   FLAG = ISLINPHASE(D) returns TRUE if the digital filter, D, is linear
%   phase. You design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   FLAG = ISLINPHASE(...,TOL) uses tolerance, TOL, to determine when two
%   numbers are close enough to be considered equal. If not specified, TOL
%   defaults to eps^(2/3).
%
%   % Example 1:
%   %   Design a Butterworth lowpass IIR filter using second order sections 
%   %   and determine if it is linear phase.
%
%   [z,p,k] = butter(6,0.7);
%   SOS = zp2sos(z,p,k);       
%   flag = islinphase(SOS)  % check if filter is linear phase
% 
%   % Example 2:
%   %   Design a highpass FIR filter and determine if it is linear phase.
%
%   b = fir1(80,0.7,'high');
%   flag = islinphase(b)    % check if filter is linear phase
% 
%   See also FVTOOL, FILTORD, ISALLPASS, ISMAXPHASE, ISMINPHASE, ISSTABLE

%   Copyright 2012-2013 The MathWorks, Inc.


narginchk(1,3);

a = 1; % Assume FIR for now

if all(size(b)>[1 1])
  % Input is a matrix, check if it is a valid SOS matrix
  if size(b,2) ~= 6
    error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
  end
    
  [b,a] = sos2tf(b); % get transfer function  

elseif ~isempty(varargin)
  a = varargin{1};
  if all(size(a)>[1 1])
    error(message('signal:signalanalysisbase:inputnotsupported'));
  end
  varargin(1) = [];
end

if isempty(varargin)
  tol = [];
else
  tol = varargin{1};
end

f = logical(signalpolyutils('islinphase',b,a,tol));

