function [T,Uinit] = tucker_als(X,R,opts)
%TUCKER_ALS Higher-order orthogonal iteration.
%
%   T = TUCKER_ALS(X,R) computes the best rank(R1,R2,..,Rn)
%   approximation of tensor X, according to the specified dimensions
%   in vector R.  The input X can be a tensor, sptensor, ktensor, or
%   ttensor.  The result returned in T is a ttensor.
%
%   T = TUCKER_ALS(X,R,OPTS) specify options:
%   OPTS.tol: Tolerance on difference in fit {1.0e-4}
%   OPTS.maxiters: Maximum number of iterations {50}
%   OPTS.dimorder: Order to loop through dimensions {1:ndims(A)}
%   OPTS.init: Initial guess [{'random'}|'eigs'|cell array]
%
%   [T,U0] = TUCKER_ALS(...) also returns the initial guess.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   T = tucker_als(X,2);        %<-- best rank(2,2,2) approximation 
%   T = tucker_als(X,[2 2 1]);  %<-- best rank(2,2,1) approximation 
%   T = tucker_als(X,2,struct('dimorder',[3 2 1]));
%   T = tucker_als(X,2,struct('dimorder',[3 2 1],'init','eigs'));
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of T
%   T = tucker_als(X,2,struct('dimorder',[3 2 1],'init',{U0}));
%
%   See also TTENSOR, TENSOR, SPTENSOR, KTENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2006, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2006) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: tucker_als.m,v 1.3 2006/12/01 19:31:39 tgkolda Exp $

% Fill in optional variable
if ~exist('opts','var')
    opts = struct;
end

% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

% Set algorithm parameters from input or by using defaults
fitchangetol = setparam(opts,'tol',1e-4);
maxiters = setparam(opts,'maxiters',50);
dimorder = setparam(opts,'dimorder',1:N);
init = setparam(opts,'init','random');

if numel(R) == 1
    R = R * ones(N,1);
end
U = cell(N,1);

%% Error checking 
% Error checking on maxiters
if maxiters < 0
    error('OPTS.maxiters must be positive');
end

% Error checking on dimorder
if ~isequal(1:N,sort(dimorder))
    error('OPTS.dimorder must include all elements from 1 to ndims(X)');
end

%% Set up and error checking on initial guess for U.
if iscell(init)
    Uinit = init;
    if numel(Uinit) ~= N
        error('OPTS.init does not have %d cells',N);
    end
    for n = dimorder(2:end);
        if ~isequal(size(Uinit{n}),[size(X,n) R(n)])
            error('OPTS.init{%d} is the wrong size',n);
        end
    end
else
    % Observe that we don't need to calculate an initial guess for the
    % first index in dimorder because that will be solved for in the first
    % inner iteration.
    if strcmp(init,'random')
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            Uinit{n} = rand(size(X,n),R(n));
        end
    elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
        % Compute an orthonormal basis for the dominant
        % Rn-dimensional left singular subspace of
        % X_(n) (1 <= n <= N).
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            fprintf('  Computing %d leading e-vectors for factor %d.\n', ...
                    R(n),n);
            Uinit{n} = nvecs(X,n,R(n));
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;

fprintf('\nAlternating Least-Squares:\n');

%% Main Loop: Iterate until convergence
for iter = 1:maxiters

    fitold = fit;

    % Iterate over all N modes of the tensor
    for n = dimorder(1:end)
	Utilde = ttm(X, U, -n, 't');

	% Maximize norm(Utilde x_n W') wrt W and
	% keeping orthonormality of W
        U{n} = nvecs(Utilde,n,R(n));
    end

    % Assemble the current approximation
    core = ttm(Utilde, U, n, 't');
    T = ttensor(core, U);

    % Compute fit
    normresidual = sqrt( normX^2 + norm(T)^2 - 2 * innerprod(X,T) );
    fit = 1 - (normresidual / normX); %fraction explained by model
    fitchange = abs(fitold - fit);

    fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);

    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end

end

%% Compute the final result
% Create the core array
core = ttm(X, U, 't');

% Assemble the resulting tensor
T = ttensor(core, U);

end

%%
function x = setparam(opts,name,default)
if isfield(opts,name);
    x = opts.(name);
else
    x = default;
end
end
