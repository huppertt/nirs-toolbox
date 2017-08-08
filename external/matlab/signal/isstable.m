function f = isstable(b,varargin)
%ISSTABLE True for stable filter
%   FLAG = ISSTABLE(B,A) returns a logical output, FLAG, equal to TRUE if
%   the filter specified by numerator coefficients B, and denominator
%   coefficients A, is stable. Input vectors B, and A define a filter with
%   transfer function:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   FLAG = ISSTABLE(SOS) returns TRUE if the filter specified using the
%   second order sections matrix, SOS, is stable. SOS is a Kx6 matrix,
%   where the number of sections, K, must be greater than or equal to 2.
%   Each row of SOS corresponds to the coefficients of a second order
%   filter. From the transfer function displayed above, the ith row of the
%   SOS matrix corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   FLAG = ISSTABLE(D) returns TRUE if the digital filter, D, is stable. 
%   You design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   % Example 1:
%   %   Create an unstable filter and verify its instability.
%
%   b = [1 2 3 4 5 5 1 2];      % numerator coefficients
%   a = [4 5 6 7 9 10 4 6];     % denominator coefficients
%   flag = isstable(b,a)        % determine if the filter is stable
%   zplane(b,a)                 % zero-pole plot for filter
%             
%   % Example 2:
%   %   Create a filter and determine its stability for different 
%   %   coefficient data types and tolerances.  
%
%   b = [1 -.5];                                % numerator coefficients
%   a = [1 -.999999999];                        % denominator coefficients
%   act_flag1 = isstable(b,a)                   % determine if its stable
%   act_flag2 = isstable(single(b),single(a))   % becomes unstable due to 
%                                               % precision
%   zplane(b,a)                                 % zero-pole plot for filter
% 
%   % Example 3:
%   %   Design a Butterworth highpass IIR filter using second order sections 
%   %   and determine its stability.
%
%   [z,p,k] = butter(6,0.7,'high');
%   SOS = zp2sos(z,p,k);    
%   flag = isstable(SOS)        % determine if the filter is stable
%   zplane(z,p)                 % zero-pole plot for filter
%
%   See also FVTOOL, FILTORD, ISALLPASS, ISLINPHASE, ISMAXPHASE, ISMINPHASE

%   Copyright 2012-2013 The MathWorks, Inc.

narginchk(1,2);

isTF = true;

if all(size(b)>[1 1])
  % Input is a matrix, check if it is a valid SOS matrix
  if size(b,2) ~= 6
    error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
  end
  if nargin > 1
    error(message('signal:signalanalysisbase:invalidNumInputs'));
  end
  isTF = false;
else
  if isempty(varargin)
    % Only numerator has been passed, so this is an FIR filter which is
    % always stable. 
    f = true;
    return;
  else
    a = varargin{1};    
    if all(size(a)>[1 1])
      error(message('signal:signalanalysisbase:inputnotsupported'));
    end
  end
end

if isTF
  if logical(signalpolyutils('isfir',b,a))
    % Section is FIR, always stable
    f = true;
  else
    f = logical(signalpolyutils('isstable',a));
  end
else
  f = true;
  for indx = 1:size(b,1)
    f = f && logical(signalpolyutils('isstable',b(indx,4:6)));
  end
end


