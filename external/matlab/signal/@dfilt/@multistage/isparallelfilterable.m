function logi = isparallelfilterable(this)
%ISPARALLELFILTERABLE   

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

logi = true;
for n = 1:length(this.Stage),
    logi = logi && isparallelfilterable(this.Stage(n));
end

% [EOF]
