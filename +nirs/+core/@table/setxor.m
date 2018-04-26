function [c,ia,ib] = setxor(a,b,varargin)
%SETXOR Find rows that occur in one or the other of two tables, but not both.
%   C = SETXOR(A,B) for tables A and B, returns the set of rows that are not in
%   the intersection of the two arrays, with repetitions removed. The rows in
%   the table C are sorted.
%
%   A and B must have the same variable names, except for order.  INTERSECT(A,B)
%   works on complete rows of A and B, considering all of their variables.  To
%   find the exclusive or with respect to a subset of those variables, use
%   column subscripting such as SETXOR(A(:,VARS),B(:,VARS)), where VARS is a
%   positive integer, a vector of positive integers, a variable name, a cell
%   array of variable names, or a logical vector.
%
%   SETXOR does not take row names into account, two rows that have the same
%   values but different names are considered equal.
% 
%   [C,IA,IB] = SETXOR(A,B) also returns index vectors IA and IB such that C is
%   a sorted combination of the values A(IA,:) and B(IB,:). If there are
%   repeated rows in A or B, then the index of the first occurrence is returned.
%
%   [C,...] = SETXOR(A,B,'stable') returns the rows in C in the same order that
%   they appear in A, then B.
%
%   [C,...] = SETXOR(A,B,'sorted') returns the rows in C in sorted order.  This
%   is the default behavior.
%
%   See also UNION, INTERSECT, SETDIFF, ISMEMBER,
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

% Calling setxor with either 'sorted' or 'stable' gives occurrence='first'
[d,ia,ib] = setxorLocal(ainds,binds,flag,'rows');

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
function [d,ia,ib] = setxorLocal(a,b,order,~)
% The main function doesn't actually need the rows themselves, since those
% are just dummy indices anyway.  It needs to know which of the two inputs
% each row of the xor "came from", so rather than returning the rows, this
% local function returns a logical indicating rows of the result that came
% from the second input (true), or from the first (false).

% Make sure a and b contain unique elements.
[uA,ia] = unique(a,'rows',order);
[uB,ib] = unique(b,'rows',order);

catuAuB = [uA;uB];                    % Sort [uA;uB] in order to find matching entries
[sortuAuB,indSortuAuB] = sortrows(catuAuB);
d = find(all(sortuAuB(1:end-1,:)==sortuAuB(2:end,:),2));    % d indicates the location of matching entries
indSortuAuB([d;d+1]) = [];                                  % Remove all matching entries - indSortuAuB only contains elements not in intersect
if strcmp(order, 'stable') % 'stable'
    indSortuAuB = sort(indSortuAuB);  % Sort the indices to get 'stable' order
end

n = size(uA,1);
d = indSortuAuB > n;           
ia = ia(indSortuAuB(~d,1),1);      % Find indices in indSortuAuB that belong to A 
ib = ib(indSortuAuB(d,1)-n,1);     % Find indices in indSortuAuB that belong to B

