function z = eq(x,y)
%EQ Equal (==) for sptensors.
%
%   A == B compares the elements of A and B for equality. The arguments can
%   be a pair of sptensors, an sptensor and a tensor, or an sptensor and a
%   scalar.  Regardless, the result is always returned as a sparse tensor.
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
% $Id: eq.m,v 1.1 2006/11/14 16:46:00 tgkolda Exp $

%% Case 1a: x is a sparse tensor and y is a scalar
if isa(x,'sptensor') && isnumeric(y) && numel(y) == 1
    if y == 0        
        z = ~x;
    else
        idx = (x.vals == y);
        z = sptensor(x.subs(idx,:),true,size(x));
    end   
    return;
end

%% Case 1b: x is a scalar and y is a sparse tensor
% Call back with the arguments reversed.
if isa(y,'sptensor') && isnumeric(x) && numel(x) == 1
    z = eq(y,x);
    return;
end

%% Case 2: Both x and y are tensors or some sort
% Check that the sizes match
if ~isequal(x.size,y.size)
    error('Size mismatch');
end

%% Case 2a: Two sparse tensors
if isa(x,'sptensor') && isa(y,'sptensor')

    % Find where their zeros intersect
    xzerosubs = setdiff(allsubs(x),x.subs,'rows');
    yzerosubs = setdiff(allsubs(y),y.subs,'rows');
    zzerosubs = intersect(xzerosubs,yzerosubs,'rows');
    
    % find where their nonzeros intersect 
    [nzsubs,ix,iy] = intersect(x.subs,y.subs,'rows');
    znzsubs = nzsubs(x.vals(ix) == y.vals(iy),:);

    % Build z
    z = sptensor([zzerosubs;znzsubs],true,x.size);
    
    return;

end

%% Case 2b: y is a dense tensor
if isa(y,'tensor')

    % Find where their zeros intersect
    yzerosubs = find(y == 0);
    zzerosubs = yzerosubs(extract(x,yzerosubs) == 0,:);

    % Find where their nonzeros intersect 
    yvals = y(x.subs);
    znzsubs = x.subs(yvals == x.vals,:);
    
    % Build z
    z = sptensor([zzerosubs;znzsubs],true,x.size);
    
    return;
    
end

%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');