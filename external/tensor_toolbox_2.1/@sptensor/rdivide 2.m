function C = rdivide(A,B)
%RDIVIDE Element-wise multiplication for sparse tensors.
%
%   RDIVIDE(A,B) denotes element-by-element division.  
% 
%   C = RDIVIDE(A,B) is called for the syntax 'A ./ B' when A is an
%   sptensor. 
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
% $Id: rdivide.m,v 1.6 2006/08/30 00:30:20 bwbader Exp $

% Scalar case
if (prod(size(B)) == 1)
    C = A;
    C.vals = C.vals/B;
    return;
end

if ~isequal(size(A),size(B))
    error('Must be two tensors of the same size');
end

switch class(B)
    case {'sptensor'}
        [tf,loc] = ismember(A.subs,B.subs,'rows');
        nzidx = find(tf);
        cvals(nzidx) = A.vals(nzidx) ./ B.vals(loc(nzidx));
        zridx = find(tf == 0);
        if ~isempty(zridx)
            warning('Divide by zero.');
        end
        cvals(zridx) = Inf;
        csubs = A.subs;
        C = sptensor(csubs, cvals, size(A));
        return;
    case {'tensor'}
        csubs = A.subs;
        cvals = A.vals ./ B(csubs); 
        C = sptensor(csubs, cvals, size(A));
        return;       
    case {'ktensor'}    
        R = numel(B.lambda);
        N = ndims(A);
        NZ = nnz(A);
        csubs = A.subs;
        avals = A.vals;       
        bvals = zeros(NZ,1);
        for r = 1:R
            tvals = B.lambda(r) * ones(NZ,1);
            for n = 1:N
                v = B{n}(:,r);
                tvals = tvals .* v(csubs(:,n));
            end
            bvals = bvals + tvals;
        end
        cvals = avals ./ bvals;
        C = sptensor(csubs, cvals, size(A));
        return;       
    otherwise
        error('Invalid second argument for sptensor/times');
end
