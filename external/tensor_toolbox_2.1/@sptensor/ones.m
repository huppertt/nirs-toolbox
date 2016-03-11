function t = ones(t)
%ONES Replace nonzero elements of sparse tensor with ones.
%
%   S = ONES(T) generates a sparse tensor with the same sparsity
%   structure as T, but with ones in the nonzero position.
%
%   See also SPTENSOR, SPONES.
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
% $Id: ones.m,v 1.8 2006/08/21 21:04:39 bwbader Exp $

t.vals = ones(size(t.vals));
