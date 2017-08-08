function s = obj2struct(this)
%OBJ2STRUCT <short description>

%   Copyright 2010 The MathWorks, Inc.

s = get(this);
if ~isempty(this.from)
    s.from = get(this.from);
end 

% [EOF]
