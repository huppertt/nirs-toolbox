function ison = parseOnOff(onoff,paramname)
% Process argument that can be true, false, 0, 1, 'on', or 'off', and
% return true if it is 'on', 1 or true.

%   Copyright 2011-2013 The MathWorks, Inc.

if     isscalar(onoff) && islogical(onoff)
    ison = onoff;
elseif isscalar(onoff) && isnumeric(onoff) && (onoff==0 || onoff==1)
    ison = logical(onoff);
elseif internal.stats.isString(onoff)
    onoff = internal.stats.getParamVal(onoff,{'on' 'off'},paramname);
    ison = strcmp(onoff,'on');
else
    m = message('stats:internal:parseOnOff:InvalidInput', paramname);
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end
