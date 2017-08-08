function st = convertScoreType(st)

%   Copyright 2014 The MathWorks, Inc.

if ~ischar(st)
    error(message('stats:classreg:learning:internal:convertScoreType:BadScoreType'));
end

allowed = {'probability' '01' 'inf' 'unknown' 'none'};
tf = strncmpi(st,allowed,length(st));
if sum(tf)~=1
    error(message('stats:classreg:learning:internal:convertScoreType:BadScoreValue',sprintf(' ''%s''',allowed{:})));
end

st = allowed{tf};
if strcmp(st,'unknown') || strcmp(st,'none')
    st = [];
end

end
