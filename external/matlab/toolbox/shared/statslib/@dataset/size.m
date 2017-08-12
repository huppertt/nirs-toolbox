function varargout = size(a,dim)
%SIZE Size of a dataset array.
%   D = SIZE(A) returns the two-element row vector D = [NOBS,NVARS] containing
%   the number of observations and variables in the dataset A.
%
%   [NOBS,NVARS] = SIZE(A) returns the number of observations and variables in
%   the dataset A as separate output variables. 
%
%   [M1,M2,M3,...,MN] = SIZE(A), for N>1, returns the sizes of the first N 
%   dimensions of the dataset A.  If the number of output arguments N does not
%   equal NDIMS(A), then for:
%
%   N > NDIMS(A), SIZE returns ones in the "extra" variables, i.e., outputs
%                 NDIMS(A)+1 through N.
%   N < NDIMS(A), MN contains the product of the sizes of dimensions N
%                 through NDIMS(A).
%
%   M = SIZE(A,DIM) returns the length of the dimension specified by the
%   scalar DIM.  For example, SIZE(A,1) returns the number of observations. If
%   DIM > NDIMS(A), M will be 1.
%
%   See also DATASET/LENGTH, DATASET/NDIMS.

%   Copyright 2006-2011 The MathWorks, Inc.


if nargin == 1
    if nargout < 2
        varargout = {[a.nobs a.nvars]};
    elseif nargout == 2
        varargout = {a.nobs a.nvars};
    else
        varargout(1:2) = {a.nobs a.nvars};
        varargout(3:nargout) = {1};
    end
else % if nargin == 2
    if nargout > 1
        error(message('stats:dataset:size:TooManyOutputs'));
    elseif ~isscalar(dim) || (dim < 1 || 2^31 < dim) || (round(dim) ~= dim)
        error(message('stats:dataset:size:InvalidDim'));
    elseif dim == 1
        varargout = {a.nobs};
    elseif dim == 2
        varargout = {a.nvars};
    else
        varargout = {1};
    end
end
