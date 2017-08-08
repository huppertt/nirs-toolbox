function f = isallpass(b,varargin)
%ISALLPASS  True for all-pass filter
%   FLAG = ISALLPASS(B,A) returns a logical output, FLAG, equal to TRUE if
%   the filter specified by numerator coefficients B, and denominator
%   coefficients A, is all-pass. Input vectors B, and A define a filter
%   with transfer function:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   FLAG = ISALLPASS(SOS) returns TRUE if the filter specified using the
%   second order sections matrix, SOS, is all-pass. SOS is a Kx6 matrix,
%   where the number of sections, K, must be greater than or equal to 2.
%   Each row of SOS corresponds to the coefficients of a second order
%   filter. From the transfer function displayed above, the ith row of the
%   SOS matrix corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   FLAG = ISALLPASS(D) returns TRUE if the digital filter, D, is all-pass.
%   You design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   FLAG = ISALLPASS(...,TOL) uses tolerance, TOL, to determine when two
%   numbers are close enough to be considered equal. If not specified, TOL
%   defaults to eps^(2/3).
%
%   % Example 1:
%   %   Create an allpass filter and verify that it is indeed allpass.
%
%   b = [1/3 1/4 1/5 1];             % numerator coefficients
%   a = fliplr(b);                   % denominator coefficients
%   flag_isallpass1 = isallpass(b,a) % determine whether filter is allpass
%   fvtool(b,a)                      % visualize the filter
% 
%   % Example 2:
%   %   Create a lattice allpass filter and verify that it is indeed allpass.
%
%   k = [1/2 1/3 1/4 1/5];          % lattice coefficients
%   [b,a] = latc2tf(k,'allpass');   % transfer function form
%   flag_isallpass = isallpass(b,a) % determine whether filter is allpass
%   fvtool(b,a)                     % visualize the filter 
% 
%   % Example 3:
%   %   Design a Butterworth highpass IIR filter using second order sections 
%   %   and verify that it is not an allpass filter.
%
%   [z,p,k] = butter(6,0.7,'high');
%   SOS = zp2sos(z,p,k);    
%   flag_isallpass = isallpass(SOS) % determine whether filter is allpass
%   fvtool(SOS)                     % visualize the filter
%
%   See also FVTOOL, FILTORD, ISLINPHASE, ISMAXPHASE, ISMINPHASE, ISSTABLE

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

f = logical(signalpolyutils('isallpass',b,a,tol));
