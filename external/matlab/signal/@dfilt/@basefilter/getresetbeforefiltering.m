function resetbeforefiltering = getresetbeforefiltering(h, dummy)
%GETRESETBEFOREFILTERING   Get the resetbeforefiltering.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if h.PersistentMemory,
    resetbeforefiltering = 'off';
else
    resetbeforefiltering = 'on';
end


% [EOF]
