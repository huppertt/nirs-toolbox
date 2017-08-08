function f = isminphase(b,varargin)
%ISMINPHASE  True for minimum phase filter
%   FLAG = ISMINPHASE(B,A) returns a logical output, FLAG, equal to TRUE if
%   the filter specified by numerator coefficients B, and denominator
%   coefficients A, is minimum phase. Input vectors B, and A define a
%   filter with transfer function:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   FLAG = ISMINPHASE(SOS) returns TRUE if the filter specified using the
%   second order sections matrix, SOS, is minimum phase. SOS is a Kx6
%   matrix, where the number of sections, K, must be greater than or equal
%   to 2. Each row of SOS corresponds to the coefficients of a second order
%   filter. From the transfer function displayed above, the ith row of the
%   SOS matrix corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   FLAG = ISMINPHASE(D) returns TRUE if the digital filter, D, is minimum
%   phase. You design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   FLAG = ISMINPHASE(...,TOL) uses tolerance, TOL, to determine when two
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
%   %   denominator coefficients, check if the filter is minimum phase for  
%   %   different tolerances.
%
%   b = single([1 1.00001]);            % numerator coefficients
%   a = single([1 .45]);                % denominator coefficients
%   min_flag1 = isminphase(b,a)         % default tolerance 
%   min_flag2 = isminphase(b,a,1e-3)    % larger tolerance
% 
%   % Example 3:
%   %   Design a Lowpass Butterworth IIR filter using second order sections
%   %   and check if it is minimum phase. 
%
%   [z,p,k] = butter(6,0.15);
%   SOS = zp2sos(z,p,k);            
%   min_flag = isminphase(SOS)  
%
%   See also FVTOOL, FILTORD, ISALLPASS, ISLINPHASE, ISMAXPHASE, ISSTABLE

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
  f = logical(signalpolyutils('isminphase',b,tol)) && ...
    logical(signalpolyutils('isstable',a));       
else
  % If SOS matrix, then check section by section  
  f = true;
  for indx = 1:size(b,1)
    f = f && logical(signalpolyutils('isminphase', b(indx, 1:3), tol)) ...
    && logical(signalpolyutils('isstable', b(indx, 4:6)));     
  end
end
