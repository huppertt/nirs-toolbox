function varargout = size(t,dim)
%SIZE Size of a table.
%   D = SIZE(T) returns the two-element row vector D = [NROWS,NVARS] containing
%   the number of rows and variables in the table T.
%
%   [NROWS,NVARS] = SIZE(T) returns the number of rows and variables in the
%   table T as separate output variables.
%
%   [M1,M2,M3,...,MN] = SIZE(T), for N>1, returns the sizes of the first N 
%   dimensions of the table T.  If the number of output arguments N does not
%   equal NDIMS(T), then for:
%
%   N > NDIMS(T), SIZE returns ones in the "extra" variables, i.e., outputs
%                 NDIMS(T)+1 through N.
%   N < NDIMS(T), MN contains the product of the sizes of dimensions N
%                 through NDIMS(T).
%
%   M = SIZE(T,DIM) returns the length of the dimension specified by the
%   scalar DIM.  For example, SIZE(T,1) returns the number of rows. If
%   DIM > NDIMS(T), M will be 1.
%
%   See also HEIGHT, WIDTH, NUMEL, NDIMS.

%   Copyright 2012-2013 The MathWorks, Inc.

if nargin == 1
    if nargout < 2
        varargout = {[t.nrows t.nvars]};
    elseif nargout == 2
        varargout = {t.nrows t.nvars};
    else
        varargout(1:2) = {t.nrows t.nvars};
        varargout(3:nargout) = {1};
    end
else % if nargin == 2
    nargoutchk(0,1);
    if ~isscalar(dim) || (dim < 1 || 2^31 < dim) || (round(dim) ~= dim)
        error(message('MATLAB:table:size:InvalidDim'));
    elseif dim == 1
        varargout = {t.nrows};
    elseif dim == 2
        varargout = {t.nvars};
    else
        varargout = {1};
    end
end
