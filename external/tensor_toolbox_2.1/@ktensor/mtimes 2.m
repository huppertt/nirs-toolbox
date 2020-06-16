function C = mtimes(A,B)
%MTIMES Implement A*B (scalar multiply) for ktensor.
%
%   C = mtimes(A,B) computes A * B where A is a Kruskal tensor and B is
%   a scalar (or vice versa). The result C is the same size as A.
%
%   See also KTENSOR.
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
% $Id: mtimes.m,v 1.4 2006/08/21 21:04:39 bwbader Exp $

% Note: We can do scalar times a tensor, but anything more complex is
% an error.

if isa(B,'numeric') && isequal(size(B),[1 1])
    C = ktensor(B * A.lambda, A.u);
elseif isa(A,'numeric') && isequal(size(A),[1 1])
    C = ktensor(A * B.lambda, B.u);
else
    error('Use mtimes(full(A),full(B)).');
end
