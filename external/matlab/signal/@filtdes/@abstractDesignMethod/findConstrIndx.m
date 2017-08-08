function indx = findConstrIndx(h,ft)
%FINDCONSTR Find the appropriate constructor.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

s = get(h,'availableTypes');

indx = 0;
true = 0;
while ~true,
    indx = indx + 1;
    true = strcmp(s(indx).tag,ft);
end

