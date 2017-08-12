function name = matchpropertyname(a,name,exact)
%MATCHPROPERTYNAME Validate a dataset array property name.

%   Copyright 2006-2011 The MathWorks, Inc.


% This matches names against this list of property names, including
% 'ObsNames'and 'VarNames', even though they are not under the 'props' field.
propertyNames = [ dataset.propsFieldNames; {'ObsNames'; 'VarNames'} ];

if ~(ischar(name) && isvector(name) && (size(name,1)==1))
    error(message('stats:dataset:matchpropertyname:InvalidPropertyName'));
end

if nargin < 3 || ~exact
    j = find(strncmp(name,propertyNames,length(name)));
else
    j = find(strcmp(name,propertyNames));
end
if isempty(j)
    error(message('stats:dataset:matchpropertyname:UnknownProperty', name));
elseif ~isscalar(j)
    error(message('stats:dataset:matchpropertyname:AmbiguousProperty', name));
end

name = propertyNames{j};
