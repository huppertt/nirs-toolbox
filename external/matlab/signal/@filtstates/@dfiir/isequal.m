function flag = isequal(h,h2)
%ISEQUAL   True if objects are numerically equal.

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.

flag = true;
flds = fieldnames(get(h));
for n = 1:length(flds),
    flag = isequal(h.(flds{n}),h2.(flds{n}));
end

% [EOF]
