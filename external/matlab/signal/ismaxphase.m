function f = ismaxphase(b,varargin)
%ISMAXPHASE  True for maximum phase filter
%   FLAG = ISMAXPHASE(B,A) returns a logical output, FLAG, equal to TRUE if
%   the filter specified by numerator coefficients B, and denominator
%   coefficients A, is maximum phase. Input vectors B, and A define a
%   filter with transfer function:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   FLAG = ISMAXPHASE(SOS) returns TRUE if the filter specified using the
%   second order sections matrix, SOS, is maximum phase. SOS is a Kx6
%   matrix, where the number of sections, K, must be greater than or equal
%   to 2. Each row of SOS corresponds to the coefficients of a second order
%   filter. From the transfer function displayed above, the ith row of the
%   SOS matrix corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   FLAG = ISMAXPHASE(D) returns TRUE if the digital filter, D, is maximum
%   phase. You design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   FLAG = ISMAXPHASE(...,TOL) uses tolerance, TOL, to determine when two
%   numbers are close enough to be considered equal. If not specified, TOL
%   defaults to eps^(2/3).
%
%   % Example 1:
%   %   Design maximum-phase and minimum-phase lattice filters and verify
%   %   their phase type.
%
%   k = [1/6 1/1.4];                % lattice coefficients
%   bmax = latc2tf(k,'max');        % convert to transfer function form
%   bmin = latc2tf(k,'min');        % convert to transfer function form
%   max_flag = ismaxphase(bmax)     % determine if bmax filter is maximum phase
%   min_flag = isminphase(bmin)     % determine if bmin filter is minimum phase
% 
%   % Example 2:
%   %   For a filter defined with a set of single precision numerator and 
%   %   denominator coefficients, check if the filter is maximum phase for  
%   %   different tolerances.
%
%   b = single([1 -0.9999]);            % numerator coefficients
%   a = single([1 0.45]);               % denominator coefficients
%   max_flag1 = ismaxphase(b,a)         % default tolerance 
%   max_flag2 = ismaxphase(b,a,1e-3)    % larger tolerance
% 
%   % Example 3:
%   %   Design a Lowpass Butterworth IIR filter using second order sections 
%   %   and check if it is maximum phase. 
%
%   [z,p,k] = butter(6,0.15);
%   SOS = zp2sos(z,p,k);            
%   max_flag = ismaxphase(SOS) 
%
%   See also FVTOOL, FILTORD, ISALLPASS, ISLINPHASE, ISMINPHASE, ISSTABLE

%   Copyright 2012-2013 The MathWorks, Inc.

narginchk(1,3);

a = 1; % Assume FIR for now
isTF = true;

if all(size(b)>[1 1])
  % Input is a matrix, check if it is a valid SOS matrix
  if size(b,2) ~= 6
    error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
  end
      
  isTF = false; % input is SOS matrix, not a transfer function

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

if isTF
  f = logical(signalpolyutils('ismaxphase',b,a,tol));
else
  % If SOS matrix, then check section by section
  f = true;  
  for indx = 1:size(b,1)
    f = f && logical(signalpolyutils('ismaxphase', b(indx,1:3), b(indx,4:6), tol));
  end  
end
