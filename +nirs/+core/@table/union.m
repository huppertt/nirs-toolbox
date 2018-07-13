function [c,ia,ib] = union(a,b,varargin)
%UNION Find rows that occur in either of two tables.
%   C = UNION(A,B) for tables A and B, returns the combined set of rows from the
%   two arrays, with repetitions removed. The rows in the table C are sorted.
%
%   A and B must have the same variable names, except for order.  UNION(A,B)
%   works on complete rows of A and B, considering all of their variables.  To
%   find the union with respect to a subset of those variables, use column
%   subscripting such as UNION(A(:,VARS),B(:,VARS)), where VARS is a positive
%   integer, a vector of positive integers, a variable name, a cell array of
%   variable names, or a logical vector.
%
%   UNION does not take row names into account, two rows that have the same
%   values but different names are considered equal.
% 
%   [C,IA,IB] = UNION(A,B) also returns index vectors IA and IB such that C is a
%   sorted combination of the values A(IA,:) and B(IB,:). If there are common
%   rows in A and B, then only the index from A is returned, in IA. If there are
%   repeated rows in A or B, then the index of the first occurrence is returned.
%
%   [C,...] = UNION(A,B,'stable') returns the rows in C in the same order that
%   they appear in A, then B.
%
%   [C,...] = UNION(A,B,'sorted') returns the rows in C in sorted order.  This
%   is the default behavior.
%
%   See also INTERSECT, SETDIFF, SETXOR,
%            ISMEMBER, UNIQUE, SORTROWS.

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

% Calling union with either 'sorted' or 'stable' gives occurrence='first'
[d,ia,ib] = unionLocal(ainds,binds,flag,'rows');

c = [subsrefParens(a,substruct('()',{ia ':'})); ...
     subsrefParens(b,substruct('()',{ib ':'}))];
% vertcat automatically merges the properties

if strcmp(flag,'sorted')
    ord(~d) = 1:length(ia);
    ord(d) = length(ia) + (1:length(ib));
    c = subsrefParens(c,substruct('()',{ord ':'}));
    % vertcat creates default row names but when we reorder for 'sorted', they
    % have the wrong row number labels.
    if isempty(a.rownames)
        if ~isempty(b.rownames)
            c.rownames(~d) = matlab.internal.table.dfltRowNames(find(~d)); %#ok<FNDSB>
        end
    elseif isempty(b.rownames)
        c.rownames(d) = matlab.internal.table.dfltRowNames(find(d)); %#ok<FNDSB>
    end
end

%-----------------------------------------------------------------------
function [d,ia,ib] = unionLocal(a,b,order,~)
% The main function doesn't actually need the rows themselves, since those
% are just dummy indices anyway.  It needs to know which of the two inputs
% each row of the union "came from", so rather than returning the rows,
% this local function returns a logical indicating rows of the result that
% came from the second input (true), or from the first (false).
[~,ndx] = unique([a;b],order,'rows');
n = size(a,1);
d = ndx > n;
ia = ndx(~d,1);
ib = ndx(d,1) - n;
