function name = matchPropertyName(name,exact)
%MATCHPROPERTYNAME Validate a table property name.

%   Copyright 2012-2013 The MathWorks, Inc.

if ~(ischar(name) && isvector(name) && (size(name,1)==1))
    error(message('MATLAB:table:InvalidPropertyName'));
end

if nargin < 2 || ~exact
    j = find(strncmpi(name,table.propertyNames,length(name)));
else
    j = find(strcmp(name,table.propertyNames));
end
if isempty(j)
    error(message('MATLAB:table:UnknownProperty', name));
elseif ~isscalar(j)
    error(message('MATLAB:table:AmbiguousProperty', name));
end

name = table.propertyNames{j};
