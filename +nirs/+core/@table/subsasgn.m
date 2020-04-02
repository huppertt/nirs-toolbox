function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignment to a table.
%   T = SUBSASGN(T,S,B) is called for the syntax T(I,J)=B, T{I,J}=B, or
%   T.VAR=B when T is a table.  S is a structure array with the fields:
%       type -- string containing '()', '{}', or '.' specifying the
%               subscript type.
%       subs -- Cell array or string containing the actual subscripts.
%
%   T(I,J) = B assigns the contents of the table B to a subset of the rows
%   and variables in the table T.  I and J are positive integers, vectors
%   of positive integers, row/variable names, cell arrays containing one or
%   more row/variable names, or logical vectors.  The assignment does not
%   use row names, variable names, or any other properties of B to modify
%   properties of T; however properties of T are extended with default
%   values if the assignment expands the number of rows or variables in T.
%   Elements of B are assigned into T by position, not by matching names.
%
%   T{I,J} = B assigns the values in the array B into elements of the table
%   T.  I and J are positive integers, or logical vectors.  Columns of B
%   are cast to the types of the target variables if necessary.  If the
%   table elements already exist, T{I,J} may also be followed by further
%   subscripting as supported by the variable.
%
%   T.VAR = B, T.(VARNAME) = B, or T.(VARINDEX) = B assigns B to a table
%   variable.  VAR is a variable name literal, VARNAME is a character
%   variable containing a variable name, or VARINDEX is a positive integer.
%   If the table variable already exists, the assignment completely
%   replaces that variable.  To assign into an element of the variable,
%   T.VAR, T.(VARNAME), or T.(VARINDEX) may be followed by further
%   subscripting as supported by the variable.  In particular,
%   T.VAR(ROWNAMES,...) = B, T.VAR{ROWNAMES,...} = B, etc. (when supported
%   by VAR) provide assignment into a table variable using row names.
%
%   T.PROPERTIES.PROPERTYNAME = P assigns to a table property.  PROPERTYNAME
%   is 'RowNames', 'VariableNames', 'Description', 'VariableDescriptions',
%   'VariableUnits', 'DimensionNames', or 'UserData'.  To assign into an element
%   of the property, T.PROPERTIES.PROPERTYNAME may also be followed by further
%   subscripting as supported by the property.
%
%   LIMITATIONS:
%  
%      You cannot assign multiple values into a table variable or property using
%      assignments such as
%         [A.CellVar{1:2}] = deal(B1,B2),
%         [A.StructVar(1:2).field] = deal(B1,B2), or
%         [A.Properties.ObsNames{1:2}] = deal(B1,B2)
%      Use multiple assignments of the form A.CellVar{1} = B1 instead.
%
%   See also TABLE, SUBSREF.

%   Copyright 2012-2013 The MathWorks, Inc.

try
    creating = isequal(t,[]);
    if creating
        t = table;
    end

    switch s(1).type
    case '()'
        t = subsasgnParens(t,s,b,creating);
    case '{}'
        %assert(~creating)
        t = subsasgnBraces(t,s,b);
    case '.'
        %assert(~creating)
        t = subsasgnDot(t,s,b);
    end
catch ME
    throwAsCaller(ME)
end