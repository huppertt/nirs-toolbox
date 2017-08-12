function [c,ia] = setdiff(a,b,vars,varargin)
%SETDIFF Set difference for dataset array observations.
%   C = SETDIFF(A,B) for dataset arrays A and B, returns the set of
%   observations that are in A but not in B, with repetitions removed.  The
%   observations in the dataset array C will be in sorted order.
% 
%   C = SETDIFF(A,B,VARS) for dataset arrays A and B, returns the set of
%   observations in A that are not present in B, considering only the
%   variables specified in VARS, with repetitions removed. The observations in
%   the dataset array C will be sorted by those variables. The values for
%   variables not specified in VARS for each observation in C are taken from
%   the corresponding observation in A. If there are multiple observations in
%   A that correspond to an observation in C, those values are taken from the
%   first occurrence.
% 
%   [C,IA] = SETDIFF(A,B,...) also returns an index vector IA such that C =
%   A(IA,:).  If there are repeated observations in A, then the index of the
%   first occurrence is returned.
%
%   [C,...] = SETDIFF(A,B,VARS,'stable') returns the observations in C in the
%   same order that they appear in A. Specify VARS as [] to use its default
%   value of all variables.
%
%   [C,...] = SETDIFF(A,B,VARS,'sorted') returns the observations in C in
%   sorted order.  This is the default behavior. Specify VARS as [] to use its
%   default value of all variables.
%
%   See also DATASET/UNION, DATASET/INTERSECT, DATASET/SETXOR, DATASET/ISMEMBER,
%            DATASET/UNIQUE, DATASET/SORTROWS.

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

% Calling setdiff with either 'sorted' or 'stable' gives occurrence='first'
[~,ia] = setdiff(ainds,binds,flag,'rows');

c = subsrefParens(a,substruct('()',{ia ':'}));
