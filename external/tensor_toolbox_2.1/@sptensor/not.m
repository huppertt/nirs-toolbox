function y = not(x)
%NOT Logical NOT (~) for sptensors.
%
%   ~X performs a logical not on the input tensor X. The result always
%   returned as a sparse tensor.
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
% $Id: not.m,v 1.1 2006/11/14 16:46:00 tgkolda Exp $

%% Then compute those indicies that are not in x
subs = setdiff(allsubs(x),x.subs,'rows');

%% Assemble final result
y = sptensor(subs,true,x.size);
