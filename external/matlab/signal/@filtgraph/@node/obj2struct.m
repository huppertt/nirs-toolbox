function s = obj2struct(this)
%OBJ2STRUCT <short description>

%   Copyright 2010 The MathWorks, Inc.

N = length(this);
for I = 1:N
    s(I) = get(this(I));
    s(I).block = obj2struct(this(I).block);
end

% [EOF]
