function defaulttw(this,band)
%DEFAULTTW   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

if band>5,
    set(this, 'TransitionWidth', .5/band);
end


% [EOF]
