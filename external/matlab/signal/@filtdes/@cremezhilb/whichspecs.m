function ws = whichspecs(h)
%WHICHSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

ws = ft_whichspecs(h);

ws(1).defval = [2400 21600];

% [EOF]
