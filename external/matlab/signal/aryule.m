function [a,e,k] = aryule( x, p)
%ARYULE   AR parameter estimation via Yule-Walker method.
%   A = ARYULE(X,ORDER) returns the coefficients of the autoregressive (AR)
%   parametric signal model estimate of X using the Yule-Walker
%   (autocorrelation) method.  This method solves the Yule-Walker equations
%   by means of the Levinson-Durbin recursion. The model has order ORDER,
%   and the output array A has ORDER+1 columns.  The coefficients along the
%   Nth row of A model the Nth column of X.  If X is a vector then A is a
%   row vector.
%
%   [A,E] = ARYULE(...) returns the final prediction error E (the variance
%   estimate of the white noise input to the AR model).
%
%   [A,E,K] = ARYULE(...) returns the vector K of reflection coefficients.
%   
%   % Example:
%   %   Estimate model order using decay of reflection coefficients.
%
%   rng default;
%   y=filter(1,[1 -0.75 0.5],0.2*randn(1024,1));
%
%   % Create AR(2) process
%   [ar_coeffs,NoiseVariance,reflect_coeffs]=aryule(y,10);
%
%   % Fit AR(10) model
%   stem(reflect_coeffs); axis([-0.05 10.5 -1 1]);
%   title('Reflection Coefficients by Lag'); xlabel('Lag');
%   ylabel('Reflection Coefficent');
%
%   See also PYULEAR, ARMCOV, ARBURG, ARCOV, LPC, PRONY.

%   Ref: S. Orfanidis, OPTIMUM SIGNAL PROCESSING, 2nd Ed.
%              Macmillan, 1988, Chapter 5
%        M. Hayes, STATISTICAL DIGITAL SIGNAL PROCESSING AND MODELING, 
%              John Wiley & Sons, 1996, Chapter 8

%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(2,2)

% Cast to enforce precision rules
p = signal.internal.sigcasttofloat(p,'double','aryule','ORDER',...
  'allownumeric');
% Checks if X is valid data
isXsingle = signal.internal.sigcheckfloattype(x,'single','aryule','X');

% perform column vector conversion before checking if 2D matrix.
if isvector(x)
  x = x(:);
end
validateattributes(x,{'numeric'},{'nonempty','finite','2d'},'arburg','X');

if size(x,1) < p
   if isvector(x)
     error(message('signal:aryule:InvalidVectorLength'));
   else
     error(message('signal:aryule:InvalidNumberOfRows'));
   end
elseif isempty(p) || ~(p == round(p))
   error(message('signal:aryule:MustBeInteger'))
end
if issparse(x)
   error(message('signal:aryule:Sparse'))
end

a = zeros(size(x,2),p+1);
e = zeros(1,size(x,2));
k = zeros(p,size(x,2));

for chan = 1:size(x,2)
  R = xcorr(x(:,chan),p,'biased');
  % LEVINSON does not support 'R' of data type 'Single'. Hence this cast.
  R = double(R);
  [a(chan,:),e(chan),k(:,chan)] = levinson(R(p+1:end),p);
end

% Cast to enforce precision rules
if isXsingle
  a = single(a);
  e = single(e);
  k = single(k);
end



