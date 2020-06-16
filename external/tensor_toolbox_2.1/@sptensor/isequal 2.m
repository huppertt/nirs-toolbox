function z = isequal(x,y)
%ISEQUAL for sptensors.
%
%   ISEQUAL(A,B) compares the sparse tensors A and B for equality.
%
%   See also SPTENSOR.
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
% $Id: isequal.m,v 1.1 2006/11/26 07:31:47 tgkolda Exp $

z = false;
if isa(x,'sptensor') && isa(y,'sptensor') && isequal(x.size,y.size) 
    z = (nnz(x-y) == 0);
end