function [c,ia,ib] = intersect(a,b,vars,varargin)
%INTERSECT Set intersection for dataset array observations.
%   C = INTERSECT(A,B) for dataset arrays A and B, returns the common set of
%   observations from the two arrays, with repetitions removed. The
%   observations in the dataset array C will be in sorted order.
% 
%   C = INTERSECT(A,B,VARS) for dataset arrays A and B, returns the set of
%   common observations from the two arrays, considering only the variables
%   specified in VARS, with repetitions removed. The observations in the
%   dataset array C will be sorted by those variables. The values for
%   variables not specified in VARS for each observation in C are taken from
%   the corresponding observation in A. If there are multiple observations in
%   A that correspond to an observation in C, those values are taken from the
%   first occurrence.
% 
%   [C,IA,IB] = INTERSECT(A,B,...) also returns index vectors IA and IB
%   such that C = A(IA,:) and C(:,VARS) = B(IB,VARS).  If there are repeated
%   observations in A or B, then the index of the first occurrence is returned.
%
%   [C,...] = INTERSECT(A,B,VARS,'stable') returns the observations in C in
%   the same order that they appear in A. Specify VARS as [] to use its
%   default value of all variables.
%
%   [C,...] = INTERSECT(A,B,VARS,'sorted') returns the observations in C in
%   sorted order.  This is the default behavior. Specify VARS as [] to use its
%   default value of all variables.
%
%   See also DATASET/UNION, DATASET/SETDIFF, DATASET/SETXOR, DATASET/ISMEMBER,
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

% Calling intersect with either 'sorted' or 'stable' gives occurrence='first'
[~,ia,ib] = intersect(ainds,binds,flag,'rows');

c = subsrefParens(a,substruct('()',{ia ':'}));
