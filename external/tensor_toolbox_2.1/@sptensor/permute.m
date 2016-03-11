function t = permute(t,order)
%PERMUTE Rearrange the dimensions of a sparse tensor.
%
%   B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they
%   are in the order specified by the vector ORDER. The result has the
%   same values of A, but the order of the subscripts needed to access
%   any particular element are rearranged as specified by ORDER.
%
%   See also SPTENSOR, PERMUTE.
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
% $Id: permute.m,v 1.4 2006/08/21 21:04:39 bwbader Exp $

% Error checking
if (ndims(order) ~= 2) || (size(order,1) ~= 1) 
    error('ORDER must be a row vector');
end
   
% Check that the permuation is valid
if ~isequal(sort(order),1:ndims(t))
    error('Invalid permutation.');
end

% Do the permutation
t.subs = t.subs(:,order);
t.size = t.size(order);
