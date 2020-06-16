function A = mtimes(A,B)
%MTIMES Scalar product for tensors.
%
%   Z = MTIMES(X,Y) multiplies the tensor X by the scalar Y (or
%   vice versa).
%
%   MTIMES(X,Y) is called for the syntax 'X * Y' when X or Y is a
%   tensor.
% 
%   Examples
%   X = tensor(rand(3,4,2));
%   W = 5 * X;
%
%   See also TENSOR, TENSOR/TTT.
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
% $Id: mtimes.m,v 1.6 2006/08/28 23:05:55 tgkolda Exp $

%% Figure out which is the tensor --- swap it to be A.
if ~isa(A,'tensor')
    tmp = B;
    B = A;
    A = tmp;
end

%% Error check on B --- it must be a scalar
if numel(B) ~= 1
    error('mtimes only supports scalar multiplication');
end

%% Do computation
A.data = B * A.data;


