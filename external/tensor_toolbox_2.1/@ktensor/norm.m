function nrm = norm(A)
%NORM Frobenius norm of a ktensor.
%
%   NORM(T) returns the Frobenius norm of a ktensor.
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
% $Id: norm.m,v 1.7 2006/08/23 16:43:24 bwbader Exp $

% Retrieve the factors of A
U = A.u;

% Compute the matrix of correlation coefficients
coefMatrix = A.lambda * A.lambda';
for i = 1:ndims(A)
  coefMatrix = coefMatrix .* (U{i}'*U{i});
end

nrm = sqrt(sum(coefMatrix(:)));

return;
