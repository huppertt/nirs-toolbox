function [B,G] = sgolay(k,F,varargin)
%SGOLAY Savitzky-Golay Filter Design.
%   B = SGOLAY(K,F) designs a Savitzky-Golay (polynomial) FIR smoothing
%   filter B.  The polynomial order, K, must be less than the frame size,
%   F, and F must be odd.  
%
%   Note that if the polynomial order K equals F-1, no smoothing
%   will occur.
%
%   SGOLAY(K,F,W) specifies a weighting vector W with length F
%   containing real, positive valued weights employed during the
%   least-squares minimization.
%
%   [B,G] = SGOLAY(...) returns the matrix G of differentiation filters.
%   Each column of G is a differentiation filter for derivatives of order
%   P-1 where P is the column index.  Given a length F signal X, an
%   estimate of the P-th order derivative of its middle value can be found
%   from:
%
%                     ^(P)
%                     X((F+1)/2) = P!*G(:,P+1)'*X
%
%   % Example:
%   %   Use sgolay to smooth a noisy sinusoid and compare both.
%
%   N = 4;                 % Order of polynomial fit
%   F = 21;                % Window length
%   [b,g] = sgolay(N,F);   % Calculate S-G coefficients
%   dx = .2;
%   xLim = 200;
%   x = 0:dx:xLim-1;
%   y = 5*sin(0.4*pi*x)+randn(size(x));  % Sinusoid with noise
%   HalfWin  = ((F+1)/2) -1;
%   SG0 = zeros(1,985); % Preallocating for speed
%   for n = (F+1)/2:996-(F+1)/2,
%   % Zero-th derivative (smoothing only)
%   SG0(n) =   dot(g(:,1), y(n - HalfWin: n + HalfWin));
%   end    
%   plot([y(1:length(SG0))', SG0'])
%   legend('Noisy Sinusoid','S-G Smoothed sinusoid')
%
%   See also SGOLAYFILT, FIR1, FIRLS, FILTER

%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(2,3);
% Cast to enforce Precision Rules
k = signal.internal.sigcasttofloat(k,'double','sgolay','K',...
  'allownumeric');
F = signal.internal.sigcasttofloat(F,'double','sgolay','F',...
  'allownumeric');

% Check if the input arguments are valid
if round(F) ~= F
  error(message('signal:sgolay:FrameMustBeInteger'))
end
if rem(F,2) ~= 1
 error(message('signal:sgolay:InvalidDimensions'))
end
if round(k) ~= k
  error(message('signal:sgolay:DegreeMustBeInteger'))
end
if k > F-1
  error(message('signal:sgolay:DegreeGeLength'))
end
if nargin < 3
   % No weighting matrix, make W an identity
   % Cast to enforce Precision Rules
   W = eye(F,'double');
else
   % Cast to enforce Precision Rules
   W = varargin{1};
   W = signal.internal.sigcasttofloat(W,'double','sgolay','W',...
     'allownumeric');
   % Check W is real.
   if ~isreal(W)
     error(message('signal:sgolay:NotReal'))
   end
   % Check for right length of W
   if length(W) ~= F
     error(message('signal:sgolay:MismatchedDimensions'))
   end
   % Check to see if all elements are positive
   if min(W) <= 0
     error(message('signal:sgolay:WVMustBePos'))
   end
   % Diagonalize the vector to form the weighting matrix
   W = diag(W);
end

% Compute the projection matrix B
s = fliplr(vander(-(F-1)/2:(F-1)/2));
S = s(:,1:k+1);   % Compute the Vandermonde matrix

[Q,R] = qr(sqrt(W)*S,0); %#ok

G = S*inv(R)*inv(R)'; % Find the matrix of differentiators

B = G*S'*W; 


% [EOF] - sgolay.m

