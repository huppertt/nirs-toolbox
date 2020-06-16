function t = full(t)
%FULL Convert a ktensor to a (dense) tensor.
%
%   T = FULL(C) converts a ktensor to a (dense) tensor.
%
%   Examples
%   X = ktensor([3; 2], rand(4,2), rand(5,2), rand(3,2));
%   Y = full(A) %<-- equivalent dense tensor
%
%   See also KTENSOR, TENSOR.
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
% $Id: full.m,v 1.10 2006/08/31 22:11:42 bwbader Exp $

sz = size(t);
data = t.lambda' * khatrirao(t.u,'r')';
t = tensor(data,sz);
