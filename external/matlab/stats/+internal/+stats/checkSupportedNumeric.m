function checkSupportedNumeric(name,x,okInteger)
%   checkSupportedNumeric(NAME,VAL) checks that the input named NAME has a
%   value VAL of a numeric type supported by code in the
%   Statistics and Machine Learning Toolbox
%
%   checkSupportedNumeric(NAME,VAL,TRUE) marks integer values as okay.

%   Copyright 2014 The MathWorks, Inc.

m = [];
if ~isreal(x)
    m = message('stats:internal:utils:NoComplexNamed',name);
elseif issparse(x)
    m = message('stats:internal:utils:NoSparseNamed',name);
elseif isobject(x)
    m = message('stats:internal:utils:NoObjectsNamed',name);
elseif (nargin<3 || ~okInteger) && ~isfloat(x)
    m = message('stats:internal:utils:FloatRequiredNamed',name);
end
if ~isempty(m)
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end
end

