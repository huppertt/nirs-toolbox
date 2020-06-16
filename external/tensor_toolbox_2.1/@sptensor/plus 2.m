function C = plus(A,B)
%PLUS Binary addition for sparse tensors. 
%
%   PLUS(A,B) adds sparse tensors A and B.  A and B must have the same
%   dimensions.  
% 
%   C = PLUS(A,B) is called for the syntax 'A + B' when A or B is a
%   tensor.
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
% $Id: plus.m,v 1.5 2006/08/21 21:04:39 bwbader Exp $

if ~isa(A,'sptensor') || ~isa(B,'sptensor') ...
	|| ~isequal(size(A),size(B))

    error('Must be two sparse tensors of the same size');

end

C = sptensor([A.subs; B.subs], [A.vals; B.vals], size(A));
