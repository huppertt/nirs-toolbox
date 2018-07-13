function [lia,locb] = ismember(a,b,varargin)
%ISMEMBER Find rows in one table that occur in another table.
%   LIA = ISMEMBER(A,B) for tables A and B returns a vector containing true in
%   elements corresponding to rows in A that are also present in B, and false
%   otherwise.
%
%   A and B must have the same variable names, except for order.  ISMEMBER(A,B)
%   works on complete rows of A and B, considering all of their variables.  To
%   find A's rows in B with respect to a subset of those variables, use column
%   subscripting such as ISMEMBER(A(:,VARS),B(:,VARS)), where VARS is a positive
%   integer, a vector of positive integers, a variable name, a cell array of
%   variable names, or a logical vector.
%
%   ISMEMBER does not take row names into account, two rows that have the same
%   values but different names are considered equal.
%
%   [LIA,LOCB] = ISMEMBER(A,B) also returns a vector LOCB containing the index
%   to the first row in B that corresponds to each row in A, or 0 if there is no
%   such row.
%
%   See also UNION, INTERSECT, SETDIFF, SETXOR,
%            UNIQUE, SORTROWS.

%   Copyright 2012-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings

if ~isa(a,'nirs.core.table') || ~isa(b,'nirs.core.table')
    error(message('MATLAB:table:setmembership:TypeMismatch'));
end

[tf,b2a] = ismember(b.varnames,a.varnames);
if ~all(tf) || (length(tf) ~= a.nvars)
    error(message('MATLAB:table:setmembership:DisjointVars'));
end
avars = 1:a.nvars;
bvars(1,b2a) = 1:b.nvars;

if nargin > 2
    narginchk(2,4);
    if ~isStrings(varargin)
        error(message('MATLAB:table:setmembership:UnknownInput2'));
    end

    % Ignore 'rows', it's always implied, but accepted anyway.  Do not accept
    % 'R2012a' or 'legacy', or 'stable' and 'sorted', or anything else.
    if any(strcmpi('legacy',varargin)) || any(strcmpi('R2012a',varargin))
        error(message('MATLAB:table:setmembership:BehaviorFlags'));
    end
    varargin(strcmpi(varargin,'rows')) = [];
    if ~isempty(varargin)
        error(message('MATLAB:table:setmembership:UnknownFlag2',varargin{1}));
    end
end

[ainds,binds] = nirs.core.table.table2midx(a,avars,b,bvars);

% Calling ismember with 'R2012a' gives occurrence='first'
[lia,locb] = ismember(ainds,binds,'rows','R2012a');
