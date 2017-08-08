function opts = getfvtooloptions(d)
%GETFVTOOLOPTIONS Returns the options sent to FVTool from the design method

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

opts = {'Analysis', 'magnitude'};

if isdb(d),
    opts = {opts{:}, 'MagnitudeDisplay', 'magnitude (db)'};
else
    opts = {opts{:}, 'MagnitudeDisplay', 'zero-phase'};
end

% [EOF]
