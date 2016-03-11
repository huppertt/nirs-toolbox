function K = fixsigns(K)
%FIXSIGNS Fix sign ambiguity of a ktensor.
%
%   K = FIXSIGNS(K) makes it so that the largest magnitude entries for
%   each vector in each factor of K are positive, provided that the
%   sign on *pairs* of vectors in a rank-1 component can be flipped.
%
%   See also KTENSOR and KTENSOR/ARRANGE.
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
% $Id: fixsigns.m,v 1.6 2006/08/25 20:36:47 bwbader Exp $

R = length(K.lambda);
for r = 1 : R
    
    for n = 1:ndims(K)
        [val(n),idx(n)] = max(abs(K.u{n}(:,r)));    
        sgn(n) = sign(K.u{n}(idx(n),r));
    end

    negidx = find(sgn == -1);
    nflip = 2 * floor(numel(negidx)/2);

    for i = 1:nflip
        n = negidx(i);
        K.u{n}(:,r) =  -K.u{n}(:,r);
    end

end
