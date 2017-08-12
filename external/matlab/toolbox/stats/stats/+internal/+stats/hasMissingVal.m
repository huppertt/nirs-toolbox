function tf = hasMissingVal(v,allowND)
% hasMissingVal Any row of this array has a missing value
%    TF = hasMissingVal(X) returns a logical array with one value per row
%    of X, indicating whether any X value in that row is missing. Missing
%    values are NaN for floating point types, "<undefined>" for categorical
%    types, empty for cell arrays of strings, or entirely blank for
%    character arrays.


%   Copyright 2013 The MathWorks, Inc.

if nargin>1
    tf = statslib.internal.hasMissingVal(v,allowND);
else
    tf = statslib.internal.hasMissingVal(v);
end

