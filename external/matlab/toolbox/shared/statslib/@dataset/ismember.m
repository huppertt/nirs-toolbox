function [lia,locb] = ismember(a,b,vars,varargin)
%ISMEMBER True for dataset array observation set members.
%   LIA = ISMEMBER(A,B) for dataset arrays A and B returns a vector containing
%   true in elements corresponding to observations in A that are also present
%   in B, and false otherwise.
%
%   LIA = ISMEMBER(A,B,VARS) returns a vector containing true in elements
%   corresponding to observations in A where the values of the variables
%   specified in VARS are the same as those for one or more observations in B,
%   and false otherwise.
%
%   [LIA,LOCB] = ISMEMBER(A,B,...) also returns a vector LOCB containing the
%   index to the first observation in B that corresponds to each observation
%   in A, or 0 if there is no such observation.
%
%   See also DATASET/UNION, DATASET/INTERSECT, DATASET/SETDIFF, DATASET/SETXOR,
%            DATASET/UNIQUE, DATASET/SORTROWS.

%   Copyright 2012 The MathWorks, Inc.


if ~isa(a,'dataset') || ~isa(b,'dataset')
    error(message('stats:dataset:setmembership:TypeMismatch'));
end

if nargin < 3 || isempty(vars)
    [tf,b2a] = ismember(b.varnames,a.varnames);
    if ~all(tf) || (length(tf) ~= a.nvars)
        error(message('stats:dataset:setmembership:DisjointVars'));
    end
    avars = 1:a.nvars;
    bvars = b2a(1:b.nvars);
else
    avars = getvarindices(a,vars,false);
    bvars = getvarindices(b,vars,false);
end

if nargin > 3
    narginchk(2,5);
    % Ignore 'rows' from the flags if present.  It's always implied, but accept it anyway.
    if any(strcmpi('legacy',varargin)) || any(strcmpi('R2012a',varargin))
        error(message('stats:dataset:setmembership:BehaviorFlags'));
    end
    varargin(strcmpi(varargin,'rows')) = [];
    if ~isempty(varargin)
        flag = varargin{1};
        if ischar(flag)
            error(message('stats:dataset:setmembership:UnknownFlag',flag));
        else
            error(message('stats:dataset:setmembership:UnknownInput'));
        end
    end
end

[ainds,binds] = dataset2idx(a,avars,b,bvars);

% Calling ismember with 'R2012a' gives occurrence='lowest'
[lia,locb] = ismember(ainds,binds,'rows','R2012a');
