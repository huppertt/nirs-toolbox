function [subs,vals] = find(t)
%FIND Find subscripts of nonzero elements in a sparse tensor.
%
%   [SUBS,VALS] = FIND(T) returns the subscripts and corresponding
%   values of the nonzero elements of T.
%
%   Note that unlike the standard MATLAB find function for an array,
%   find does not return linear indices. Instead, it returns an M x N
%   array where M is the number of nonzero values and N = ndims(T).
%   Thus, I(k,:) specifies the subscript of value V(k).
%
%   See also SPTENSOR, FIND.
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
% $Id: find.m,v 1.4 2006/08/21 21:04:39 bwbader Exp $

subs = t.subs;
vals = t.vals;
