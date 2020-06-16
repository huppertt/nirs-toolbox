function V = mttkrp(X,U,n)
%MTTKRP Matricized tensor times Khatri-Rao product for tensor.
%
%   V = MTTKRP(X,U,n) efficiently calculates the matrix product of the
%   n-mode matricization of X with the Khatri-Rao product of all
%   entries in U, a cell array of matrices, except the nth.  How to
%   most efficiently do this computation depends on the type of tensor
%   involved.
%
%   See also TENSOR, TENMAT, KHATRIRAO
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
% $Id: mttkrp.m,v 1.3 2006/08/31 22:11:42 bwbader Exp $

N = ndims(X);
%Xn = tenmat(X,n);
%Xn = double(Xn);
Xn = permute(X,[n 1:n-1,n+1:N]);
Xn = reshape(Xn.data, size(X,n), prod(size(X))/size(X,n));
Z = khatrirao(U{[1:n-1,n+1:N]},'r');
V = Xn*Z;
