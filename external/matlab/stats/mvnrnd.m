function [r,T] = mvnrnd(mu,sigma,cases,T)
%MVNRND Random vectors from the multivariate normal distribution.
%   R = MVNRND(MU,SIGMA) returns an N-by-D matrix R of random vectors
%   chosen from the multivariate normal distribution with mean vector MU,
%   and covariance matrix SIGMA.  MU is an N-by-D matrix, and MVNRND
%   generates each row of R using the corresponding row of MU.  SIGMA is a
%   D-by-D symmetric positive semi-definite matrix, or a D-by-D-by-N array.
%   If SIGMA is an array, MVNRND generates each row of R using the
%   corresponding page of SIGMA, i.e., MVNRND computes R(I,:) using MU(I,:)
%   and SIGMA(:,:,I).  If the covariance matrix is diagonal, containing
%   variances along the diagonal and zero covariances off the diagonal,
%   SIGMA may also be specified as a 1-by-D matrix or a 1-by-D-by-N array,
%   containing just the diagonal. If MU is a 1-by-D vector, MVNRND
%   replicates it to match the trailing dimension of SIGMA.
%
%   R = MVNRND(MU,SIGMA,N) returns a N-by-D matrix R of random vectors
%   chosen from the multivariate normal distribution with 1-by-D mean
%   vector MU, and D-by-D covariance matrix SIGMA.
%
%   Example:
%
%      mu = [1 -1]; Sigma = [.9 .4; .4 .3];
%      r = mvnrnd(mu, Sigma, 500);
%      plot(r(:,1),r(:,2),'.');
%
%   See also MVTRND, MVNPDF, MVNCDF, NORMRND.

%   R = MVNRND(MU,SIGMA,N,T) supplies the Cholesky factor T of
%   SIGMA, so that SIGMA(:,:,J) == T(:,:,J)'*T(:,:,J) if SIGMA is a 3D array or SIGMA
%   == T'*T if SIGMA is a matrix.  No error checking is done on T.
%
%   [R,T] = MVNRND(...) returns the Cholesky factor T, so it can be
%   re-used to make later calls more efficient, although there are greater
%   efficiency gains when SIGMA can be specified as a diagonal instead.

%   Copyright 1993-2011 The MathWorks, Inc.


if nargin < 2 || isempty(mu) || isempty(sigma)
    error(message('stats:mvnrnd:TooFewInputs'));
elseif ndims(mu) > 2
    error(message('stats:mvnrnd:BadMu'));
elseif ndims(sigma) > 3
    error(message('stats:mvnrnd:BadSigma'));
end

[n,d] = size(mu);
sz = size(sigma);
if sz(1)==1 && sz(2)>1
    % Just the diagonal of Sigma has been passed in.
    sz(1) = sz(2);
    sigmaIsDiag = true;
else
    sigmaIsDiag = false;
end


% Special case: if mu is a column vector, then use sigma to try
% to interpret it as a row vector.
if d == 1 && sz(1) == n
    mu = mu';
    [n,d] = size(mu);
end

% Get size of data.
if nargin < 3 || isempty(cases)
    nocases = true; % cases not supplied
else
    nocases = false; % cases was supplied
    if n == cases
        % mu is ok
    elseif n == 1 % mu is a single row, make cases copies
        n = cases;
        mu = repmat(mu,n,1);
    else
        error(message('stats:mvnrnd:InputSizeMismatchMu'));
    end
end

% Single covariance matrix
if ndims(sigma) == 2
    % Make sure sigma is the right size
    
    if sz(1) ~= sz(2)
        error(message('stats:mvnrnd:BadCovariance2DSize'));
    elseif ~isequal(sz, [d d])
        error(message('stats:mvnrnd:InputSizeMismatch2DSigma'));
    end
    
    if nargin > 3
        % sigma has already been factored, so use it.
        r = randn(n,size(T,1)) * T + mu;
    elseif sigmaIsDiag
        % Just the diagonal of sigma has been specified.
        if any(sigma<0)
            error(message('stats:mvnrnd:BadDiagSigma'));
        end
        t = sqrt(sigma);
        if nargout>1
            T = diag(t);
        end
        r = bsxfun(@times,randn(n,d),t) + mu;
    else
        % Factor sigma using a function that will perform a Cholesky-like
        % factorization as long as the sigma matrix is positive
        % semi-definite (can have perfect correlation). Cholesky requires a
        % positive definite matrix.  sigma == T'*T
        [T,err] = cholcov(sigma);
        if err ~= 0
            error(message('stats:mvnrnd:BadCovariance2DSymPos'));
        end
        r = randn(n,size(T,1)) * T + mu;
    end
    
% Multiple covariance matrices
elseif ndims(sigma) == 3
    % mu is a single row and cases not given, rep mu out to match sigma
    if n == 1 && nocases % already know size(sigma,3) > 1
        n = sz(3);
        mu = repmat(mu,n,1);
    end
    
    % Make sure sigma is the right size
    if sz(1) ~= sz(2) % Sigma is 3-D
        error(message('stats:mvnrnd:BadCovariance3DSize'));
    elseif (sz(1) ~= d) || (sz(2) ~= d) % Sigma is 3-D
        error(message('stats:mvnrnd:InputSizeMismatch3DSigma'));
    elseif sz(3) ~= n
        error(message('stats:mvnrnd:InputSizeMismatchSigmaDimension'));
    end
  
    if nargin < 4
        if nargout > 1
            T = zeros(sz,class(sigma));
        end
        if sigmaIsDiag
            sigma = reshape(sigma,sz(2),sz(3))';
            if any(any(sigma<0))
                error(message('stats:mvnrnd:BadDiagSigma'));
            end
            R = sqrt(sigma);
            r = bsxfun(@times,randn(n,d),R) + mu;
            if nargout > 1
                for i=1:n
                    T(:,:,i) = diag(R(i,:));
                end
            end
        else
            r = zeros(n,d,superiorfloat(mu,sigma));
            for i = 1:n
                [R,err] = cholcov(sigma(:,:,i));
                if err ~= 0
                    error(message('stats:mvnrnd:BadCovariance3DSymPos'));
                end
                Rrows = size(R,1);
                r(i,:) = randn(1,Rrows) * R + mu(i,:);
                if nargout > 1
                    T(1:Rrows,:,i) = R;
                end
            end
        end
    else
        % T specified
        r = zeros(n,d,superiorfloat(mu,sigma));
        for i = 1:n
            r(i,:) = randn(1,d) * T(:,:,i) + mu(i,:);
        end
    end
end
