function siz = size(a,idx)
%SIZE Return size of sptenmat.
%  
%   D = SIZE(T) returns the size of the tensor.  
%
%   I = size(T,DIM) returns the sizes of the dimensions specified by DIM.
%
%   See also SPTENMAT.
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
% $Id: size.m,v 1.5 2006/08/21 21:04:39 bwbader Exp $

if isempty(a.tsize)
    siz = [];
    return;
end

m = prod(a.tsize(a.rdims));
n = prod(a.tsize(a.cdims));
siz = [m n];

if  exist('idx','var')
    siz = siz(idx);
end
