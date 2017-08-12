function a = subsasgn(a,s,b)
%SUBSASGN Subscripted assignment to a dataset array.
%   A = SUBSASGN(A,S,B) is called for the syntax A(I,J)=B, A{I,J}=B, or
%   A.VAR=B when A is a dataset array.  S is a structure array with the
%   fields:
%       type -- string containing '()', '{}', or '.' specifying the
%               subscript type.
%       subs -- Cell array or string containing the actual subscripts.
%
%   A(I,J) = B assigns the contents of the dataset array B to a subset of the
%   observations and variables in the dataset array A.  I and J are positive
%   integers, vectors of positive integers, observation/variable names, cell
%   arrays containing one or more observation/variable names, or logical
%   vectors.  The assignment does not use observation names, variable names,
%   or any other properties of B to modify properties of A; however properties
%   of A are extended with default values if the assignment expands the number
%   of observations or variables in A. Elements of B are assigned into A by
%   position, not by matching names.
%
%   A{I,J} = B assigns the value B into an element of the dataset array A.  I
%   and J are positive integers, or logical vectors.  Cell indexing cannot
%   assign into multiple dataset elements, that is, the subscripts I and J
%   must each refer to only a single observation or variable.  B is cast to
%   the type of the target variable if necessary.  If the dataset element
%   already exists, A{I,J} may also be followed by further subscripting as
%   supported by the variable.
%
%   For dataset variables that are cell arrays, assignments such as
%   A{1,'CellVar'} = B assign into the contents of the target dataset element
%   in the same way that {}-indexing of an ordinary cell array does.
%
%   For dataset variables that are N-D arrays, i.e., each observation is a
%   matrix or array, assignments such as A{1,'ArrayVar'} = B assigns into the
%   second and following dimensions of the target dataset element, i.e., the
%   assignment adds a leading singleton dimension to B to account for the
%   observation dimension of the dataset variable.
%
%   A.VAR = B or A.(VARNAME) = B assigns B to a dataset variable.  VAR is a
%   variable name literal, or VARNAME is a character variable containing a
%   variable name.  If the dataset variable already exists, the assignment
%   completely replaces that variable.  To assign into an element of the
%   variable, A.VAR or A.(VARNAME) may be followed by further subscripting as
%   supported by the variable.  In particular, A.VAR(OBSNAMES,...) = B and
%   A.VAR{OBSNAMES,...} = B (when supported by VAR) provide assignment into a
%   dataset variable using observation names.
%
%   A.PROPERTIES.PROPERTYNAME = P assigns to a dataset property.  PROPERTYNAME
%   is 'ObsNames', 'VarNames', 'Description', 'VarDescription', 'Units',
%   'DimNames', or 'UserData'.  To assign into an element of the property,
%   A.PROPERTIES.PROPERTYNAME may also be followed by further subscripting as
%   supported by the property.
%
%
%   LIMITATIONS:
%
%      You cannot assign multiple values into dataset variables or properties
%      using assignments such as [A.CellVar{1:2}] = B,
%      [A.StructVar(1:2).field] = B, or [A.Properties.ObsNames{1:2}] = B.  Use
%      multiple assignments of the form A.CellVar{1} = B instead.
%
%      Similarly, if a dataset variable is a cell array with multiple columns
%      or is an N-D cell array, then the contents of that variable for a
%      single observation consists of multiple cells, and you cannot assign to
%      all of them using the syntax A{1,'CellVar'} = B.  Use multiple
%      assignments of the form [A.CellVar{1,1}] = B instead.
%
%   See also DATASET, DATASET/SUBSREF, DATASET/SET.

%   Copyright 2006-2012 The MathWorks, Inc.


creating = isequal(a,[]);
if creating
    a = dataset;
end

switch s(1).type
case '()'
    a = subsasgnParens(a,s,b,creating);
case '{}'
    a = subsasgnBraces(a,s,b);
case '.'
    a = subsasgnDot(a,s,b);
end
