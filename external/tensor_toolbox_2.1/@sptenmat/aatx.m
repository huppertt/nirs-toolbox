function z = aatx(a,x)
%AATX Implicitly compute A * A' * x for sptenmat.
%
%   Z = AATX(A,X) takes a sptenmat object A and computes A * A' *
%   X. This is done without converting A to a standard MATLAB sparse
%   matrix.
%
%   This function is likely most useful as an argument to a routine
%   such as EIGS.
%
%   See also SPTENMAT, SPTENSOR/EIGS.
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
% $Id: aatx.m,v 1.5 2006/08/21 21:04:39 bwbader Exp $

subs = a.subs;
s1 = subs(:,1);
s2 = subs(:,2);
m = size(a,1);
n = size(a,2);
vals = a.vals;

v1 = x(s1);
v1 = vals .* v1;
y = accumarray(s2, v1, [n 1]);

v2 = y(s2);
v2 = vals .* v2;
z = accumarray(s1, v2, [m 1]);

