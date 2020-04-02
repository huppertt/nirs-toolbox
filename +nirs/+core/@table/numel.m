function n = numel(t,varargin)
%NUMEL Number of elements in a table.
%   N = NUMEL(T) returns the number of elements in the table T, equivalent to
%   PROD(SIZE(T)).  Note that variables in a table may themselves have multiple
%   columns.  NUMEL(T) does not account for that.
%
%   N = NUMEL(T, INDEX1, INDEX2, ...) returns the number of subscripted elements
%   in T(INDEX1, INDEX2, ...), equivalent to PROD(SIZE(T(INDEX1, INDEX2, ...))).
%
%   See also SIZE, HEIGHT, WIDTH.

%   Copyright 2012-2013 The MathWorks, Inc.

switch nargin
case 1
    n = t.nrows * t.nvars;

otherwise
    % Return the total number of elements in the subscript expression.  Don't do
    % any checks to see if the subscripts actually exist, just count them up.
    % subsref/subsasgn will error later on if the subscripts refer to something
    % that's not there.

    if numel(varargin) ~= t.ndims
        error(message('MATLAB:table:NDSubscript'));
    end

    rowIndices = varargin{1};
    if ischar(rowIndices)
        if strcmp(rowIndices,':') % already checked ischar
            nrows = t.nrows;
        else
            nrows = 1;
        end
    elseif islogical(rowIndices)
        nrows = nnz(rowIndices);
    elseif isnumeric(rowIndices) || iscellstr(rowIndices)
        nrows = numel(rowIndices);
    end

    varIndices = varargin{2};
    if ischar(varIndices)
        if strcmp(varIndices,':') % already checked ischar
            nvars = t.nvars;
        else
            nvars = 1;
        end
    elseif islogical(varIndices)
        nvars = nnz(varIndices);
    elseif isnumeric(varIndices) || iscellstr(varIndices)
        nvars = numel(varIndices);
    end
    n = nrows*nvars;
end
