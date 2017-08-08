function bool = thisisquantizable(Hd)
%THISISQUANTIZABLE Returns true if the dfilt object can be quantized

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

if length(Hd.Stage) > 1 & ~isa(Hd.Stage(1), 'dfilt.multistage'),
    bool = isa(Hd.Stage, class(Hd.Stage(1)));
else
    bool = false;
end

% [EOF]
