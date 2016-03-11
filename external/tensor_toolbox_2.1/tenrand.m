function X = tenrand(sz)
%TENRAND Uniformly distributed pseudo-random tensor.
%
%   X = TENRAND(SZ) forms a tensor of size SZ with pseudo-random
%   values drawn from a uniform distribution on the unit interval.
%
%   TENRAND(SZ) is equivalent to TENSOR(RAND(SZ(1),SZ(2),...),SZ).
%
%   See also TENSOR, RAND.
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
% $Id: tenrand.m,v 1.2 2006/08/21 21:04:39 bwbader Exp $

data = rand([sz 1 1]);
X = tensor(data,sz);
