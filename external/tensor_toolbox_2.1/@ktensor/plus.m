function C = plus(A,B)
%PLUS Binary addition for ktensor.
%
%   C = PLUS(A,B) adds two ktensors of the same size, and the
%   result is a ktensor of the same size.
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
% $Id: plus.m,v 1.7 2006/08/24 23:59:07 bwbader Exp $

if (isa(A,'ktensor') && isa(B,'ktensor'))    

    if ~isequal(size(A),size(B))
	error('Tensor size mismatch.')
    end

    lambda = [A.lambda; B.lambda];    
    M = ndims(A);
    u = cell(M,1);
    for m = 1 : M
        u{m} = [A.u{m} B.u{m}];
    end 
    C = ktensor(lambda, u);
    return;
end

error('Use plus(full(A),full(B)).');
