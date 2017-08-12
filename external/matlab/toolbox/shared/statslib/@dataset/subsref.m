function [varargout] = subsref(a,s)
%SUBSREF Subscripted reference for a dataset array.
%   B = SUBSREF(A,S) is called for the syntax A(I,J), A{I,J}, or A.VAR
%   when A is a dataset array.  S is a structure array with the fields:
%       type -- string containing '()', '{}', or '.' specifying the
%               subscript type.
%       subs -- Cell array or string containing the actual subscripts.
%
%   B = A(I,J) returns a dataset array that contains a subset of the
%   observations and variables in the dataset array A.  I and J are positive
%   integers, vectors of positive integers, observation/variable names, cell
%   arrays containing one or more observation/variable names, or logical
%   vectors.  B contains the same property values as A, subsetted for
%   observations or variables where appropriate.
%
%   B = A{I,J} returns an element of a dataset variable.  I and J are positive
%   integers, or logical vectors.  Cell indexing cannot return multiple
%   dataset elements, that is, the subscripts I and J must each refer to only
%   a single observation or variable.  A{I,J} may also be followed by further
%   subscripting as supported by the variable.
%
%   For dataset variables that are cell arrays, expressions such as A{1,'CellVar'}
%   return the contents of the referenced dataset element in the same way that
%   {}-indexing on an ordinary cell array does.  If the dataset variable is a
%   single column of cells, the contents of a single cell is returned.  If the
%   dataset variable has multiple columns or is N-D, multiple outputs containing
%   the contents of multiple cells are returned.
%
%   For dataset variables that are N-D arrays, i.e., each observation is a
%   matrix or an array, expressions such as A{1,'ArrayVar'} return
%   A.ArrayVar(1,:,...) with the leading singleton dimension squeezed out.
%
%   B = A.VAR or A.(VARNAME) returns a dataset variable.  VAR is a variable
%   name literal, or VARNAME is a character variable containing a variable
%   name.  A.VAR or A.(VARNAME) may also be followed by further subscripting as
%   supported by the variable.  In particular, A.VAR(OBSNAMES,...) and
%   A.VAR{OBSNAMES,...} (when supported by VAR) provide subscripting into a
%   dataset variable using observation names.
%
%   P = A.PROPERTIES.PROPERTYNAME returns a dataset property.  PROPERTYNAME is
%   'ObsNames', 'VarNames', 'Description', 'VarDescription', 'Units', 'DimNames',
%   or 'UserData'.  A.PROPERTIES.PROPERTYNAME may also be followed by further
%   subscripting as supported by the property.
%
%   LIMITATIONS:
%
%      Subscripting expressions such as A.CellVar{1:2}, A.StructVar(1:2).field,
%      or A.Properties.ObsNames{1:2} are valid, but result in SUBSREF
%      returning multiple outputs in the form of a comma-separated list.  If
%      you explicitly assign to output arguments on the LHS of an assignment,
%      for example, [cellval1,cellval2] = A.CellVar{1:2}, those variables will
%      receive the corresponding values. However, if there are no output
%      arguments, only the first output in the comma-separated list is
%      returned.
%
%      Similarly, if a dataset variable is a cell array with multiple columns
%      or is an N-D cell array, then subscripting expressions such as
%      A{1,'CellVar'} result in SUBSREF returning the contents of multiple
%      cells.  You should explicitly assign to output arguments on the LHS of
%      an assignment, for example, [cellval1,cellval2] = A{1,'CellVar'}.
%
%   See also DATASET, DATASET/SUBSASGN, DATASET/GET.

%   Copyright 2006-2012 The MathWorks, Inc.


switch s(1).type
case '()'
    [varargout{1:nargout}] = subsrefParens(a,s);
case '{}'
    [varargout{1:nargout}] = subsrefBraces(a,s);
case '.'
    [varargout{1:nargout}] = subsrefDot(a,s);
end
