function Y = squeeze(X)
%SQUEEZE Remove singleton dimensions from a sparse tensor.
%
%   Y = SQUEEZE(X) returns a sparse tensor Y with the same elements as
%   X but with all the singleton dimensions removed.  A singleton
%   is a dimension such that size(X,dim)==1.  
%
%   If X has *only* singleton dimensions, then Y is a scalar.
%
%   Examples
%   squeeze( sptensor(tensor(rand(2,1,3))) ) %<-- returns a 2-by-3 tensor
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
% $Id: squeeze.m,v 1.3 2006/08/21 21:04:39 bwbader Exp $

if ~any(X.size == 1)
  % No singleton dimensions to squeeze
  Y = X;
else
  idx = find(X.size > 1)
  if numel(idx) == 0
    % Scalar case - only singleton dimensions
    Y = X.vals;
  else
    siz = X.size(idx);
    Y = sptensor(X.subs(:,idx), X.vals, siz);
  end
end

return;
