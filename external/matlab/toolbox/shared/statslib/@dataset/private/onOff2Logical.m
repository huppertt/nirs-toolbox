function paramVal = onOff2Logical(paramVal,paramName)
%ONOFF2LOGICAL Convert 'on'/'off' to logical.

%   Copyright 2011 The MathWorks, Inc.


if ischar(paramVal)
    if strcmpi(paramVal,'on')
        paramVal = true;
    elseif strcmpi(paramVal,'off')
        paramVal = false;
    else
        error(message('stats:dataset:dataset:InvalidOnOffVal', paramName));
    end
elseif isscalar(paramVal) && (islogical(paramVal) || isnumeric(paramVal))
    % leave it alone
else
%     error(message('stats:dataset:tdfread:InvalidReadVarNames'));
    error(message('stats:dataset:dataset:InvalidOnOffVal', paramName));
end
end % function


