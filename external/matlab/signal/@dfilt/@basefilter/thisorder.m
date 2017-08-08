function n = thisorder(this)
%THISORDER   Dispatch and recall.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

n = [];
Hd = dispatch(this);
for indx = 1:length(Hd)
    n = [n thisorder(Hd(indx))];
end

% [EOF]
