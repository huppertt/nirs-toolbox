function varargout = arburg( x, p)
%ARBURG   AR parameter estimation via Burg method.
%   A = ARBURG(X,ORDER) returns the coefficients of the autoregressive (AR)
%   parametric signal model estimate of X using Burg's method. The model
%   has order ORDER, and the output array A has ORDER+1 columns.  The
%   coefficients along the Nth row of A model the Nth column of X.  If X is
%   a vector then A is a row vector.
%
%   [A,E] = ARBURG(...) returns the final prediction error E (the variance
%   estimate of the white noise input to the AR model).
%
%   [A,E,K] = ARBURG(...) returns the reflection coefficients (parcor
%   coefficients) in each column of K.
%
%   % Example:
%   %   Estimate input noise variance for AR(4) model.
%   A=[1 -2.7607 3.8106 -2.6535 0.9238]; 
%
%   % Generate noise standard deviations
%   % Seed random number generator for reproducible results
%   rng default;
%   noise_stdz=rand(1,50)+0.5;
%
%   % Generate column vectors that have corresponding standard deviation
%   x = bsxfun(@times,noise_stdz,randn(1024,50));
%
%   % filter each column using the AR model.
%   y = filter(1,A,x);
%
%   % Compute the estimated coefficients and deviations for each column
%   [ar_coeffs,NoiseVariance]=arburg(y,4);
%
%   %Display the mean value of each estimated polynomial coefficient
%   estimatedA = mean(ar_coeffs)
%
%   %Compare actual vs. estimated variances
%   plot(noise_stdz.^2,NoiseVariance,'k*');
%   xlabel('Input Noise Variance');
%   ylabel('Estimated Noise Variance');
%
%   See also PBURG, ARMCOV, ARCOV, ARYULE, LPC, PRONY.

%   Ref: S. Kay, MODERN SPECTRAL ESTIMATION,
%              Prentice-Hall, 1988, Chapter 7
%        S. Orfanidis, OPTIMUM SIGNAL PROCESSING, 2nd Ed.
%              Macmillan, 1988, Chapter 5

%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(2,2)

% Checks if X is valid data
signal.internal.sigcheckfloattype(x,'','arburg','X');

% convert to column vector before 2d check is performed
if isvector(x)
  x = x(:);
end

validateattributes(x,{'numeric'},{'nonempty','finite','2d'},'arburg','X');
validateattributes(p,{'numeric'},{'positive','integer','scalar'},'arburg','ORDER');
% Cast to enforce precision rules
p = double(p);

if issparse(x),
   error(message('signal:arburg:Sparse'))
end
if size(x,1) < p+1
   if isvector(x)
      error(message('signal:arburg:VectorTooSmall', p + 1));
   else 
      error(message('signal:arburg:MatrixTooSmall', p + 1));
   end
end
N  = size(x,1);
Nchan = size(x,2);

% Initialization
ef = x;
eb = x;
% Data type of 'a' should be the same as 'x' to enforce precision rules
a = ones(1,Nchan, class(x)); %#ok<*ZEROLIKE>

% Initial error
E = dot(x,x)./N;

% Preallocate 'k' for speed.
% Data type of 'k' should be the same as 'x' to enforce precision rules
k = zeros(p, Nchan, class(x));

for m=1:p
   % Calculate the next order reflection (parcor) coefficient
   efp = ef(2:end,:);
   ebp = eb(1:end-1,:);
   num = -2.*dot(ebp,efp);
   den = dot(efp,efp)+dot(ebp,ebp);
   
   k(m,:) = num ./ den;
   
   % Update the forward and backward prediction errors
   ef = efp + bsxfun(@times,k(m,:),ebp);
   eb = ebp + bsxfun(@times,conj(k(m,:)),efp);
   
   % Update the AR coeff.
   a=[a;zeros(1,Nchan)] + bsxfun(@times,k(m,:),[zeros(1,Nchan);conj(flipud(a))]);
   
   % Update the prediction error
   E = bsxfun(@times,1 - bsxfun(@times,conj(k(m,:)),k(m,:)),E);
end

a = a.'; % By convention all polynomials are row vectors
varargout{1} = a;
if nargout >= 2
    varargout{2} = E(end,:);
end
if nargout >= 3
    varargout{3} = k;
end
