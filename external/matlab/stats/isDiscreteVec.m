function tf = isDiscreteVec(v)
% isDiscreteVec Variable is most naturally treated as a discrete column vector.
%    TF = isDiscreteVec(X) returns TRUE if X is a character array, or if
%    X is a column vector that is categorical, a cell array of strings, or
%    logical.


%   Copyright 2011 The MathWorks, Inc.

% Require a discrete variable to be a single column
if isa(v,'categorical') || iscellstr(v) || islogical(v)
    tf = isvector(v) && (size(v,2) == 1);
elseif ischar(v)
    tf = ismatrix(v);
else
    tf = false;
end
