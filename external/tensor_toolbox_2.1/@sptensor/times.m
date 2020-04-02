function C = times(A,B)
%TIMES Element-wise multiplication for sparse tensors.
%
%   TIMES(A,B) denotes element-by-element multiplication.  
% 
%   C = TIMES(A,B) is called for the syntax 'A .* B' when A or B is a
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
% $Id: times.m,v 1.7 2006/08/23 00:23:20 tgkolda Exp $

if ~isequal(size(A),size(B))
    error('Must be two tensors of the same size');
end

switch class(B)
    case {'sptensor'}
        [csubs,ia,ib] = intersect(A.subs,B.subs,'rows');
        cvals = A.vals(ia) .* B.vals(ib);
        C = sptensor(csubs, cvals, size(A));
        return;
    case {'tensor'}
        csubs = A.subs;
        cvals = A.vals .* B(csubs); 
        C = sptensor(csubs, cvals, size(A));
        return;       
    case {'ktensor'}    
        csubs = A.subs;
        cvals = zeros(size(A.vals));       
        R = numel(B.lambda);
        N = ndims(A);
        for r = 1:R
            tvals = B.lambda(r) * A.vals;
            for n = 1:N
                v = B{n}(:,r);
                tvals = tvals .* v(csubs(:,n));
            end
            cvals = cvals + tvals;
        end
        C = sptensor(csubs, cvals, size(A));
        return;       
    otherwise
        error('Invalid second argument for sptensor/times');
end