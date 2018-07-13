function [c,ia,ib] = intersect(a,b,varargin)
%INTERSECT Find rows common to two tables.
%   C = INTERSECT(A,B) for tables A and B, returns the common set of rows from
%   the two tables, with repetitions removed.  The rows in the table C are in
%   sorted order.
%
%   A and B must have the same variable names, except for order.  INTERSECT(A,B)
%   works on complete rows of A and B, considering all of their variables.  To
%   find the intersection with respect to a subset of those variables, use
%   column subscripting such as INTERSECT(A(:,VARS),B(:,VARS)), where VARS is a
%   positive integer, a vector of positive integers, a variable name, a cell
%   array of variable names, or a logical vector.
%
%   INTERSECT does not take row names into account, two rows that have the same
%   values but different names are considered equal.
% 
%   [C,IA,IB] = INTERSECT(A,B) also returns index vectors IA and IB such that C
%   = A(IA,:) and C = B(IB,:).  If there are repeated rows in A or B, then the
%   index of the first occurrence is returned.
%
%   [C,...] = INTERSECT(A,B,'stable') returns the rows in C in the same order
%   that they appear in A.
%
%   [C,...] = INTERSECT(A,B,'sorted') returns the rows in C in sorted order.
%   This is the default behavior.
%
%   See also UNION, SETDIFF, SETXOR, ISMEMBER,
%            UNIQUE, SORTROWS.

%   Copyright 2012-2014 The MathWorks, Inc.

if ~isa(a,'table') || ~isa(b,'table')
    error(message('MATLAB:table:setmembership:TypeMismatch'));
end

[tf,b2a] = ismember(b.varnames,a.varnames);
if ~all(tf) || (length(tf) ~= a.nvars)
    error(message('MATLAB:table:setmembership:DisjointVars'));
end
avars = 1:a.nvars;
bvars(1,b2a) = 1:b.nvars;

if nargin < 3
    flag = 'sorted';
else
    narginchk(2,5); % high=5, to let setmembershipFlagChecks sort flags out
    flag = setmembershipFlagChecks(varargin);
end

[ainds,binds] = table.table2midx(a,avars,b,bvars);

% Calling intersect with either 'sorted' or 'stable' gives occurrence='first'
[~,ia,ib] = intersect(ainds,binds,flag,'rows');

c = subsrefParens(a,substruct('()',{ia ':'}));

% Use b's property values where a's were empty.
if isempty(a.rownames) && ~isempty(b.rownames)
    c.rownames = b.rownames(ib);
end
c.props = mergeProps(a.props,b.props);
