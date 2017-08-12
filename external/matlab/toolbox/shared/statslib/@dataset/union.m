function [c,ia,ib] = union(a,b,vars,varargin)
%UNION Set union for dataset array observations.
%   C = UNION(A,B) for dataset arrays A and B, returns the combined set of
%   observations from the two arrays, with repetitions removed. The
%   observations in the dataset array C are sorted.
% 
%   C = UNION(A,B,VARS) returns the combined set of observations from the two
%   arrays, with repetitions of unique combinations of the variables specified
%   in VARS removed. The observations in the dataset array C are sorted by
%   those variables. The values for variables not specified in VARS for each
%   observation in C are taken from the corresponding observation in A or B,
%   or from A if there are common observations in both A and B. If there are
%   multiple observations in A or B that correspond to an observation in C,
%   those values are taken from the first occurrence.
% 
%   [C,IA,IB] = UNION(A,B,...) also returns index vectors IA and IB such that
%   C is a sorted combination of the values A(IA,:) and B(IB,:). If there are
%   common observations in A and B, then only the index from A is returned, in
%   IA. If there are repeated observations in A or B, then the index of the
%   first occurrence is returned.
%
%   [C,...] = UNION(A,B,VARS,'stable') returns the observations in C in the
%   same order that they appear in A, then B. Specify VARS as [] to use its
%   default value of all variables.
%
%   [C,...] = UNION(A,B,VARS,'sorted') returns the observations in C in sorted
%   order.  This is the default behavior. Specify VARS as [] to use its
%   default value of all variables.
%
%   See also DATASET/INTERSECT, DATASET/SETDIFF, DATASET/SETXOR,
%            DATASET/ISMEMBER, DATASET/UNIQUE, DATASET/SORTROWS.

%   Copyright 2012 The MathWorks, Inc.


if ~isa(a,'dataset') || ~isa(b,'dataset')
    error(message('stats:dataset:setmembership:TypeMismatch'));
end

[tf,b2a] = ismember(b.varnames,a.varnames);
if ~all(tf) || (length(tf) ~= a.nvars)
    error(message('stats:dataset:setmembership:DisjointVars'));
end

if nargin < 3 || isempty(vars)
    avars = 1:a.nvars;
    bvars = b2a(1:b.nvars);
else
    avars = getvarindices(a,vars,false);
    bvars = getvarindices(b,vars,false);
end

if nargin < 4
    flag = 'sorted';
else
    narginchk(2,5);
    flag = setmembershipFlagChecks(varargin);
end

[ainds,binds] = dataset2idx(a,avars,b,bvars);

% Calling union with either 'sorted' or 'stable' gives occurrence='first'
[~,ia,ib] = union(ainds,binds,flag,'rows');

c = [subsrefParens(a,substruct('()',{ia ':'})); ...
     subsrefParens(b,substruct('()',{ib ':'}))];
