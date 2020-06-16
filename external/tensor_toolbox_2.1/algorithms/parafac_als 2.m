function [P,Uinit] = parafac_als(X,R,opts)
%PARAFAC_ALS Compute a PARAFAC decomposition of any type of tensor.
%
%   P = PARAFAC_ALS(X,R) computes an estimate of the best rank-R
%   PARAFAC model of a tensor X using an alternating least-squares
%   algorithm.  The input X can be a tensor, sptensor, ktensor, or
%   ttensor. The result P is a ktensor.
%
%   P = PARAFAC_ALS(X,R,OPTS) specify options:
%   OPTS.tol: Tolerance on difference in fit {1.0e-4}
%   OPTS.maxiters: Maximum number of iterations {50}
%   OPTS.dimorder: Order to loop through dimensions {1:ndims(A)}
%   OPTS.init: Initial guess [{'random'}|'nvecs'|cell array]
%
%   [P,U0] = PARAFAC_ALS(...) also returns the initial guess.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   P = parafac_als(X,2);
%   P = parafac_als(X,2,struct('dimorder',[3 2 1]));
%   P = parafac_als(X,2,struct('dimorder',[3 2 1],'init','nvecs'));
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
%   P = parafac_als(X,2,struct('dimorder',[3 2 1],'init',{U0}));
%
%   See also KTENSOR, TENSOR, SPTENSOR, TTENSOR.
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
% $Id: parafac_als.m,v 1.11 2006/12/01 19:31:39 tgkolda Exp $

%% Fill in optional variable
if ~exist('opts','var')
    opts = struct;
end

%% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Set algorithm parameters from input or by using defaults
fitchangetol = setparam(opts,'tol',1e-4);
maxiters = setparam(opts,'maxiters',50);
dimorder = setparam(opts,'dimorder',1:N);
init = setparam(opts,'init','random');

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
        if ~isequal(size(Uinit{n}),[size(X,n) R])
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
            Uinit{n} = rand(size(X,n),R);
        end
    elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            fprintf('  Computing %d leading e-vectors for factor %d.\n',R,n);
            Uinit{n} = nvecs(X,n,R);
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

        % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
        Unew = mttkrp(X,U,n);

        % Compute the matrix of coefficients for linear system
        Y = ones(R,R);
        for i = [1:n-1,n+1:N]
            Y = Y .* (U{i}'*U{i});
        end
        Unew = (Y \ Unew')';

        % Normalize each vector to prevent singularities in coefmatrix
        if iter == 1
            lambda = sqrt(sum(Unew.^2,1))'; %2-norm
        else
            lambda = max( max(Unew,[],1), 1 )'; %max-norm
        end
        Unew = Unew * spdiags(1./lambda,0,R,R);
        U{n} = Unew;
    end

    P = ktensor(lambda,U);
    normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
    fit = 1 - (normresidual / normX); %fraction explained by model
    fitchange = abs(fitold - fit);

    fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);

    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end

end

%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = arrange(P);
% Fix the signs
P = fixsigns(P);

end

%% 
function x = setparam(opts,name,default)
if isfield(opts,name);
    x = opts.(name);
else
    x = default;
end
end
