function C = mtimes(A,B)
%MTIMES Sparse tensor times a scalar.
%
%   MTIMES(A,B) denotes multiplication.  
% 
%   C = TIMES(A,B) is called for the syntax 'A .* B' when A or B is a
%   tensor.
%
%   See also SPTENSOR, SPTENSOR/TIMES, SPTENSOR/TTT
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
% $Id: mtimes.m,v 1.7 2006/08/21 21:04:39 bwbader Exp $

if isa(A,'sptensor') && ~isa(B,'sptensor') && ...
	(numel(B) == 1)

    C = sptensor(A.subs, A.vals * B, size(A));
    return;
end

if isa(B,'sptensor') && ~isa(A,'sptensor') && ...
	(numel(A) == 1)

    C = sptensor(B.subs, B.vals * A, size(B));
    return;
end

error('Mtimes only supports a sparse tensor times a scalar');

