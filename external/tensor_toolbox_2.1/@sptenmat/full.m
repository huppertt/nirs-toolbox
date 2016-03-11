function B = full(A)
%FULL Convert a sptenmat to a (dense) tenmat.
%
%   B = FULL(A) converts a sptenmat A to a (dense) tenmat B.
%
%   See also SPTENMAT, TENMAT.
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
% $Id: full.m,v 1.1 2006/08/22 17:12:05 bwbader Exp $

% Extract the order and size of A
siz = size(A);

% Create a dense zero tensor B that is the same size as A
B = tenmat(zeros([siz,1,1]), A.rdims, A.cdims, A.tsize);

% Extract the linear indices of entries in A
idx = tt_sub2ind(siz,A.subs);

% Copy the values of A into B using linear indices
B(idx) = A.vals;
