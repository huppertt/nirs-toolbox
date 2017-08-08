function [kr,R0]=poly2rc(a,efinal)
%POLY2RC  Convert prediction polynomial to reflection coefficients (step-down). 
%   K = POLY2RC(A) returns the reflection coefficients, K, based on the 
%   prediction polynomial, A.
%
%   If A(1) is not equal to 1, POLY2RC normalizes the prediction
%   polynomial by A(1).
%
%   [K,R0] = POLY2RC(A,Efinal) returns the zero lag autocorrelation, R0,
%   based on the final prediction error, Efinal. If Efinal is not
%   specified, then the default is Efinal=0.
%
%   % Example:
%   %   Convert the following prediction filter polynomial to reflection 
%   %   coefficients:
%   %   a = [1.0000   0.6149   0.9899   0.0000   0.0031  -0.0082];
%
%   a = [1.0000   0.6149   0.9899   0.0000   0.0031  -0.0082];
%   efinal = 0.2;               % Final prediction error
%   [k,r0] = poly2rc(a,efinal)  % Reflection coefficients
%
%   See also RC2POLY, POLY2AC, AC2POLY, RC2AC, AC2RC and TF2LATC.

%   References: S. Kay, Modern Spectral Estimation,
%               Prentice Hall, N.J., 1987, Chapter 6.
%
%   Copyright 1988-2013 The MathWorks, Inc.

if (size(a,1) > 1) && (size(a,2) > 1)
    error(message('signal:poly2rc:inputnotsupported'));
end

if nargin == 1 | isempty(efinal) %#ok
    efinal = 0;
end

% Cast to enforce Precision Rules
if any([signal.internal.sigcheckfloattype(a,'single','poly2rc','A') ...
    signal.internal.sigcheckfloattype(efinal,'single','poly2rc','Efinal')])
  a = single(a);
  efinal = single(efinal);
end

if length(a) <= 1,
  % K is length of A minus one so make empty if A is a scalar or empty.
  % Cast to enforce Precision Rules
  if isa(a,'single')
    kr = single([]);
  else
    kr = [];
  end
  R0 = efinal;
  return
end


if a(1) == 0,
    error(message('signal:poly2rc:SignalErr'));
end

a = a(:)./a(1);    % Convert to column vector and normalize by a(1)

p = length(a)-1;   % The leading one does not count
e = zeros(p,1);
kr = zeros(p,1,class(a)); %#ok<*ZEROLIKE>

e(p) = efinal;

kr(p) = a(end);

for k = p-1:-1:1,
  [a,e(k)] = levdown(a,e(k+1));
  kr(k) = a(end);
end 

% Compute R0 only if asked for because it can cause divide by zero warnings
if nargout >= 2,    
  % R0 is the zero order prediction error when the prediction error filter,
  % A(z) = 1.
  R0 = e(1)./(1-abs(kr(1))^2);
end

