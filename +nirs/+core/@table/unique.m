function [c,ia,ic] = unique(a,varargin)
%UNIQUE Find unique rows in a table.
%   C = UNIQUE(A) returns a table C that contains only the sorted unique rows in
%   A. A must contain only variables whose class has a UNIQUE method. These
%   include variables that are numeric, character, logical, categorical, or cell
%   arrays of strings. For a variable with multiple columns, its class's UNIQUE
%   method must support the 'rows' flag.
%
%   UNIQUE(A) works on complete rows of A, considering all of its variables.  To
%   find unique rows with respect to a subset of those variables, use column
%   subscripting such as UNIQUE(A(:,VARS)), where VARS is a positive integer, a
%   vector of positive integers, a variable name, a cell array of variable
%   names, or a logical vector.
%
%   UNIQUE does not take row names into account, two rows that have the same
%   values but different names are considered equal.
%
%   [C,IA,IC] = UNIQUE(A) also returns column index vectors IA and IC such that
%   C = A(IA,:) and A = C(IC,:).
%
%   [C,IA,IC] = UNIQUE(A,OCCURRENCE) specifies which index is returned in IA in
%   the case of repeated rows in A. The default value is OCCURRENCE='first', which
%   returns the index of the first occurrence of each repeated row in A, while
%   OCCURRENCE='last' returns the index of the last occurrence of each
%   repeated row in A.
%
%   [C,IA,IC] = UNIQUE(A,'stable') returns the rows of C in the same order that
%   they appear in A, while [C,IA,IC] = UNIQUE(A,'sorted') returns the rows of C
%   in sorted order.
%
%   See also UNION, INTERSECT, SETDIFF, SETXOR,
%            ISMEMBER, SORTROWS.

%   Copyright 2012-2013 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings

% Remove 'rows' from the flags if present.  It's always implied, but accept it anyway.
varargin(strcmpi('rows',varargin)) = [];

vars = 1:a.nvars;
ainds = nirs.core.table.table2midx(a,vars);

% Do not accept 'R2012a' or 'legacy'.
if nargin > 1
    if ~isStrings(varargin)
        error(message('MATLAB:table:setmembership:UnknownInput2'));
    elseif any(strcmpi('legacy',varargin)) || any(strcmpi('R2012a',varargin))
        error(message('MATLAB:table:setmembership:BehaviorFlags'));
    end
end
[~,ia,ic] = unique(ainds,'rows',varargin{:});
c = subsrefParens(a,substruct('()',{ia ':'}));
