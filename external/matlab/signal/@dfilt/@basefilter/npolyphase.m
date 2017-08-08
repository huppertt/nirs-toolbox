function n = npolyphase(this)
%NPOLYPHASE   Return the number of polyphases for this filter.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

n = [];
for indx = 1:length(this)
    n = [n thisnpolyphase(this(indx))];
end

% [EOF]
