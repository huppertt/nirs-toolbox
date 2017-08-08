function regrCheck(Y,Yfit,W)

%   Copyright 2010 The MathWorks, Inc.



if length(Y)~=length(Yfit)
    error(message('stats:classreg:learning:internal:regrCheck:LengthYandYfitMismatch'));
end
if ~isvector(Y) || ~isvector(Yfit)
    error(message('stats:classreg:learning:internal:regrCheck:BadDims'));
end

% Check types of Yfit and Y
if ~isnumeric(Y)
    error(message('stats:classreg:learning:internal:regrCheck:BadY'));
end
if ~isnumeric(Yfit)
    error(message('stats:classreg:learning:internal:regrCheck:BadYfit'));
end

% Check weights
if ~isfloat(W) || ~isvector(W) || length(W)~=length(Y) || any(W<0)
    error(message('stats:classreg:learning:internal:regrCheck:BadWeights', length( Y )));
end
end
