function idx = getfreqrange(this)
%GETFREQRANGE   Return range of object and the index to range options.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

if ishalfnyqinterval(this),    
    idx = 1;     % 0 to pi       or    0 to Fs/2

else
    idx = 2;     % 0 to 2pi      or    0 to Fs

    if this.centerdc,
        idx = 3; % -pi to pi     or    -Fs/2 to Fs/2
    end
end

% [EOF]
