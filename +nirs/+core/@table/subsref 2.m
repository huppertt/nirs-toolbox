function [varargout] = subsref(t,s)
%SUBSREF Subscripted reference for a table.
%   B = SUBSREF(T,S) is called for the syntax T(I,J), T{I,J}, or T.VAR
%   when T is a table.  S is a structure array with the fields:
%       type -- string containing '()', '{}', or '.' specifying the
%               subscript type.
%       subs -- Cell array or string containing the actual subscripts.
%
%   B = T(I,J) returns a table B that contains a subset of the rows and
%   variables in the table T.  I and J are positive integers, vectors of
%   positive integers, row/variable names, cell arrays containing one or
%   more row/variable names, or logical vectors.  B contains the same
%   property values as T, subsetted for rows or variables where
%   appropriate.
%
%   B = T{I,J} returns an array B created as the horizontal concatenation
%   of the table variables specified by J, containing only those rows
%   specified by I.  SUBSREF throws an error if the types of the variables
%   are not compatible for concatenation.  I and J are positive integers,
%   vectors of positive integers, row/variable names, cell arrays
%   containing one or more row/variable names, or logical vectors.  T{I,J}
%   may also be followed by further subscripting as supported by the
%   variable.
%
%   B = T.VAR, T.(VARNAME), or T.(VARINDEX) returns a table variable.  VAR
%   is a variable name literal, VARNAME is a character variable containing
%   a variable name, or VARINDEX is a positive integer.  T.VAR, T.(VARNAME),
%   or T.(VARINDEX) may also be followed by further subscripting as supported
%   by the variable.  In particular, T.VAR(ROWNAMES,...), T.VAR{ROWNAMES,...},
%   etc. (when supported by VAR) provide subscripting into a table variable
%   using row names.
%
%   P = T.PROPERTIES.PROPERTYNAME returns a table property.  PROPERTYNAME is
%   'RowNames', 'VariableNames', 'Description', 'VariableDescriptions',
%   'VariableUnits', 'DimensionNames', or 'UserData'.  T.PROPERTIES.PROPERTYNAME
%   may also be followed by further subscripting as supported by the property.
%
%   LIMITATIONS:
%
%      Subscripting expressions such as
%         T.CellVar{1:2},
%         T.StructVar(1:2).field, or
%         T.Properties.RowNames{1:2}
%      are valid, but result in SUBSREF returning multiple outputs in the form
%      of a comma-separated list.  If you explicitly assign to output arguments
%      on the LHS of an assignment, for example,
%         [cellval1,cellval2] = T.CellVar{1:2},
%      those variables will receive the corresponding values.  However, if there
%      are no output arguments, only the first output in the comma-separated
%      list is returned.
%
%   See also TABLE, SUBSASGN.

%   Copyright 2012-2013 The MathWorks, Inc.
try
    switch s(1).type
    case '()'
        [varargout{1:nargout}] = subsrefParens(t,s);
    case '{}'
        [varargout{1:nargout}] = subsrefBraces(t,s);
    case '.'
        [varargout{1:nargout}] = subsrefDot(t,s);
    end
catch ME
    throwAsCaller(ME)
end
