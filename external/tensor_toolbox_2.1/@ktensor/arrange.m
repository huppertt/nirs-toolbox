function X = arrange(X,N)
%ARRANGE Arranges the rank-1 terms of a ktensor.
%
%   ARRANGE(X) normalizes the columns of each matrix, absorbing the
%   excess weight into lambda and then sorts everything so that the
%   lambda values are in decreasing order.
%
%   ARRANGE(X,N) absorbs the weights into the Nth factor instead of
%   lambda.
%
%   See also KTENSOR.
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
% $Id: arrange.m,v 1.8 2006/08/23 00:20:04 tgkolda Exp $

% Ensure that matrices are normalized
for r = 1:length(X.lambda)
    for n = 1:ndims(X)
        tmp = norm(X.u{n}(:,r));
        X.u{n}(:,r) = X.u{n}(:,r) / tmp;
        X.lambda(r) = X.lambda(r) * tmp;
    end
end

% Sort
[X.lambda, idx] = sort(X.lambda, 1, 'descend');
for i = 1 : ndims(X)
    X.u{i} = X.u{i}(:,idx);
end

% Absorb the weight into one factor, if requested
if exist('N','var')
    r = length(X.lambda);
    X.u{N} = X.u{N} * spdiags(X.lambda,0,r,r);
    X.lambda = ones(size(X.lambda));
end

return;
