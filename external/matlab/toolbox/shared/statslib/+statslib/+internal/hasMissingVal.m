function tf = hasMissingVal(v,allowND)
% hasMissingVal Any row of this array has a missing value
%    TF = hasMissingVal(X) returns a logical array with one value per row
%    of X, indicating whether any X value in that row is missing. Missing
%    values are NaN for floating point types, "<undefined>" for categorical
%    types, empty for cell arrays of strings, or entirely blank for
%    character arrays.


%   Copyright 2011 The MathWorks, Inc.

if nargin > 1 && allowND
    v = v(:,:);
elseif ~ismatrix(v)
    m = message('stats:internal:hasMissingVal:NDData');
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end

if isfloat(v)
    tf = any(isnan(v),2);
elseif isa(v,'categorical')
    tf = any(isundefined(v),2);
elseif iscellstr(v)
    tf = any(cellfun('isempty',v),2);
elseif ischar(v)
    tf = all(v == ' ',2);
else
    tf = false(size(v,1),1);
end
