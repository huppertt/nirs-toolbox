function b = isspecify(h)
%ISSPECIFY Returns true if OrderMode is set to 'specify'

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdynpropenab(h, 'orderMode'),
    
    orderMode = get(h,'orderMode');
    orderModeOpts = set(h,'orderMode');
    
    switch orderMode,
        case orderModeOpts{1}, %'specify'
            b = true;
        otherwise
            b = false;
    end
else
    b = true;
end

% [EOF]
