function C = minus(A,B)
%MINUS Binary subtraction for ktensor.  
%
%   C = MINUS(A,B) computes C = A - B.  A and B must both be ktensors
%   and have the same size, and the result is another ktensor of the
%   same size.
%
%   C = MINUS(A,B) is called for the syntax 'A - B' when A or B is a
%   ktensor.
%
%   See also KTENSOR, SIZE, ISEQUAL.
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
% $Id: minus.m,v 1.7 2006/08/23 16:43:24 bwbader Exp $

if (isa(A,'ktensor') && isa(B,'ktensor'))    

    if ~isequal(size(A),size(B))
	error('Tensor size mismatch.')
    end

    lambda = [A.lambda; -B.lambda];    
    M = ndims(A);
    u = cell(M,1);
    for m = 1 : M
        u{m} = [A.u{m} B.u{m}];
    end 
    C = ktensor(lambda,u);
    return;

end

error('Use minus(full(A),full(B)).');



