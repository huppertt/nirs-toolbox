function n = ndims(t)
%NDIMS Return the number of dimensions of a tensor.
%
%   NDIMS(X) returns the number of dimensions of tensor X.
%
%   Examples
%   A = rand(4,3,1); ndims(A) %<-- Returns 2
%   X = tensor(A); ndims(X) %<-- Returns 2
%   X = tensor(A,[4 3 1]); ndims(X) %<-- Returns 3
%
%   See also TENSOR
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
% $Id: ndims.m,v 1.4 2006/08/21 21:04:40 bwbader Exp $

n = numel(t.size);
