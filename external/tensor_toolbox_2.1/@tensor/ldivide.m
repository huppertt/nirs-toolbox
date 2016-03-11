function Z = ldivide(X,Y)
%LDIVIDE Left array divide for tensor.
%
%   X.\Y denotes element-by-element division.  X and Y
%   must have the same dimensions unless one is a scalar.
%   A scalar can be divided with anything.
% 
%   LDIVIDE(X,Y) is called for the syntax 'X .\ Y' when X or Y is
%   a tensor.
%
%   See also TENSOR.
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
% $Id: ldivide.m,v 1.5 2006/08/21 21:04:40 bwbader Exp $

Z = tenfun(@ldivide,X,Y);
