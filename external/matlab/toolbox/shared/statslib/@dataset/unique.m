function [c,ia,ic] = unique(a,vars,varargin)
%UNIQUE Unique observations in a dataset array.
%   C = UNIQUE(A) returns a dataset array C that contains only the sorted
%   unique observations in A. A must contain only variables whose class has a
%   UNIQUE method. These include variables that are numeric, character,
%   logical, categorical, or cell arrays of strings. For a variable with
%   multiple columns, its class's UNIQUE method must support the 'rows' flag.
%
%   [C,IA,IC] = UNIQUE(A) also returns index vectors IA and IC such that
%   C = A(IA,:) and A = C(IC,:).
%
%   C = UNIQUE(A,VARS) returns a dataset array C that contains only one
%   observation from A for each unique combination of values for the variables
%   specified in VARS. Observations in C are sorted by the variables specified
%   in VARS. VARS is a positive integer, a vector of positive integers, a
%   variable name, a cell array containing one or more variable names, or a
%   logical vector. C includes all variables from A. The values in C for the
%   variables not specified in VARS are taken from the last occurrence among
%   observations in A for each unique combination of values for the variables
%   specified in VARS.
%
%   [C,IA,IC] = UNIQUE(A,VARS) also returns index vectors IA and IC such that
%   C = A(IA,:) and A(:,VARS) = C(IC,VARS).
%
%   [C,IA,IC] = UNIQUE(A,VARS,OCCURRENCE) specifies which index is returned in
%   IA in the case of repeated observations in A. The default value is
%   OCCURENCE='last', which returns the index of the last occurrence of each
%   repeated observation in A, while OCCURRENCE='first' returns the index of
%   the first occurrence of each repeated observation in A. The values in C
%   for variables not specified in VARS are taken from the observations
%   A(IA,:). Specify VARS as [] to use its default value of all variables.
%
%   [C,IA,IC] = UNIQUE(A,VARS,'stable') returns the observations of C in the
%   same order that they appear in A, while [C,IA,IC] = UNIQUE(A,VARS,'sorted')
%   returns the observations of C in sorted order. IA and IC are column
%   vectors. If there are repeated observations in A, then IA returns the
%   index of the first occurrence of each repeated observation. Specify VARS
%   as [] to use its default value of all variables.
% 
%   In a future release, the behavior of the following syntaxes will change
%   including:
%     -	Default occurrence of indices will switch from last to first
%     -	IA and IC will always be column index vectors
% 
%   In order to see what impact those changes will have on your code, use:
% 
%        [C,IA,IC] = UNIQUE(A,VARS,'R2012a')
%        [C,IA,IC] = UNIQUE(A,VARS,OCCURRENCE,'R2012a')
% 
%   If the changes in behavior adversely affect your code, you may preserve
%   the current behavior with:
% 
%        [C,IA,IC] = UNIQUE(A,VARS,'legacy') 
%        [C,IA,IC] = UNIQUE(A,VARS,OCCURRENCE,'legacy')
%
%   See also DATASET/UNION, DATASET/INTERSECT, DATASET/SETDIFF, DATASET/SETXOR,
%            DATASET/ISMEMBER, DATASET/SORTROWS.

%   Copyright 2006-2012 The MathWorks, Inc.


if nargin < 2 || isempty(vars)
    vars = 1:a.nvars;
else
    vars = getvarindices(a,vars,false);
end

% Remove 'rows' from the flags if present.  It's always implied, but accept it anyway.
varargin(strcmpi('rows',varargin)) = [];

ainds = dataset2idx(a,vars);

[~,ia,ic] = unique(ainds,'rows',varargin{:});
c = subsrefParens(a,substruct('()',{ia ':'}));
