function f = ispolyphase(this)
%ISPOLYPHASE   Returns true if the filter is polyphase.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Don't use BASE_IS because we DON'T want to call dispatch

for indx = 1:length(this)
    f(indx) = thisispolyphase(this(indx));
end

% Make sure that we return logicals.
f = logical(f);

% [EOF]
