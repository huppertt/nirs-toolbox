function B = full(A)
%FULL Convert a sparse tensor to a (dense) tensor.
%
%   B = FULL(A) converts a sptensor A to a (dense) tensor B.
%
%   See also SPTENSOR, TENSOR.
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
% $Id: full.m,v 1.9 2006/11/26 07:31:26 tgkolda Exp $

% Extract the order and size of A
siz = size(A);

% Create a dense zero tensor B that is the same size as A
B = tensor(zeros([siz,1,1]),siz);

if isempty(A.subs)
    return;
end

% Extract the linear indices of entries in A
idx = tt_sub2ind(siz,A.subs);

% Copy the values of A into B using linear indices
B(idx) = A.vals;

