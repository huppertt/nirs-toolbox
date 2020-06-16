function z = ge(x,y)
%GE Greater than or equal for sptensors.
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
% $Id: ge.m,v 1.2 2006/11/21 00:37:38 tgkolda Exp $

%% Case 1: One argument is a scalar
if isa(x,'sptensor') && isnumeric(y) && numel(y) == 1
    subs1 = x.subs((x.vals >= y),:);
    if y <= 0
        subs2 = setdiff(allsubs(x),x.subs,'rows');
    else
        subs2 = [];
    end
    z = sptensor([subs1;subs2],true,size(x));
    return;
end

%% Case 1: One argument is a scalar
% Call back with the arguments reversed.
if isa(y,'sptensor') && isnumeric(x) && numel(x) == 1
    z = le(y,x);
    return;
end

%% Case 2: Both x and y are tensors or some sort
% Check that the sizes match
if ~isequal(x.size,y.size)
    error('Size mismatch');
end

%% Case 2a: Two sparse tensors
if isa(x,'sptensor') && isa(y,'sptensor')
    z = le(y,x);
    return;
end

%% Case 2b: y is a dense tensor
if isa(y,'tensor')

    % x zero 
    subs1 = find(y <= 0);
    subs1 = setdiff(subs1,x.subs,'rows');
    
    % x nonzero
    subs2 = x.subs(x.vals >= y(x.subs,'extract'),:);
    
    % assemble
    z = sptensor([subs1;subs2],true,size(x));
    
    return;
    
end

%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');