function n = numel(a,varargin)
%NUMEL Number of elements in a dataset array.
%   N = NUMEL(A) returns 1.  To find the number of elements, N, in the dataset
%   array A, use PROD(SIZE(A)) or NUMEL(A,':',':').
%
%   N = NUMEL(A, VARARGIN) returns the number of subscripted elements, N, in
%   A(index1, index2, ..., indexN), where VARARGIN is a cell array whose
%   elements are index1, index2, ... indexN.
%
%   See also DATASET/SIZE, DATASET/LENGTH.

%   Copyright 2006-2011 The MathWorks, Inc.


switch nargin
case 1
    % Return a 1 so that dataset.varname can return one output.
    n = 1;
%     n = a.nobs * a.nvars;

otherwise
    % Return the total number of elements in the subscript expression.
    
    if numel(varargin) ~= a.ndims
        error(message('stats:dataset:numel:NDSubscript'));
    end

    obsIndices = varargin{1};
    if ischar(obsIndices)
        if strcmp(obsIndices,':') % already checked ischar
            nrows = a.nobs;
        else
            nrows = 1;
        end
    elseif isnumeric(obsIndices) || islogical(obsIndices) || iscellstr(obsIndices)
        nrows = numel(obsIndices);
    end

    varIndices = varargin{2};
    if ischar(varIndices)
        if strcmp(varIndices,':') % already checked ischar
            ncols = a.nvars;
        else
            ncols = 1;
        end
    elseif isnumeric(varIndices) || islogical(varIndices) || iscellstr(varIndices)
        ncols = numel(varIndices);
    end
    n = nrows*ncols;
end
