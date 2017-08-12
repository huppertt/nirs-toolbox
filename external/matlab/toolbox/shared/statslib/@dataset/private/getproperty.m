function [varargout] = getproperty(a,name,createIfEmpty)
%GETPROPERTY Get a dataset array property.

%   Copyright 2006-2011 The MathWorks, Inc.


if nargin < 3, createIfEmpty = false; end

% We may be given a name (when called from get), or a subscript expression
% that starts with a '.name' subscript (when called from subsref).  Get the
% name and validate it in any case.
if isstruct(name)
    s = name;
    if s(1).type == '.'
        name = s(1).subs;
    else
        error(message('stats:dataset:getproperty:InvalidSubscript'));
    end
    haveSubscript = true;
else
    haveSubscript = false;
end
% Allow partial match for property names if this is via the get method;
% require exact match if it is via subsref
name = matchpropertyname(a,name,haveSubscript);

% Get the property out of the dataset.  Some properties need special handling
% when empty:  create either a non-empty default version or a "canonical" 0x0
% cell array (subscripting can sometimes make them 1x0 or 0x1), depending on
% what the caller asks for.
switch name
case 'ObsNames'
    p = a.obsnames;
    if isempty(p)
        if createIfEmpty
            p = dfltobsnames(1:a.nobs);
        else
            p = {}; % force 0x0
        end
    end
case 'VarNames'
    p = a.varnames;
    % varnames are "always there", so leave them 1x0 when empty
case 'DimNames'
    p = a.props.DimNames;
case {'VarDescription' 'Units'}
    p = a.props.(name);
    if isempty(p)
        if createIfEmpty
            p = repmat({''},1,a.nvars);
        else
            p = {}; % force 0x0
        end
    end
case 'Description'
    p = a.props.Description;
case 'UserData'
    p = a.props.UserData;
end

if haveSubscript && ~isscalar(s)
    % If there's cascaded subscripting into the property, let the property's
    % subsasgn handle the reference.  This may return a comma-separated list,
    % so ask for and assign to as many outputs as we're given.  If there's no
    % LHS to the original expression (nargout==0), this only assigns one
    % output and drops everything else in the CSL.
    [varargout{1:nargout}] = subsref(p,s(2:end));
else
    % If there's no cascaded subscripting, only ever assign the property itself.
    varargout{1} = p;
end
