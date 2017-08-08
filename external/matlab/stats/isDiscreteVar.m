function tf = isDiscreteVar(v)
% isDiscreteVar Variable is most naturally treated as discrete.
%    TF = isDiscreteVar(X) returns TRUE if X is categorical, a cell array
%    of strings, logical, or a character array.


%   Copyright 2011-2014 The MathWorks, Inc.

% For most purposes a variable must be a vector. However, a char variable
% is a matrix with each row treated as a category name.
tf =    (isvector(v) && (isa(v,'categorical') || iscellstr(v) || islogical(v))) ...
     || (ismatrix(v) && ischar(v));

