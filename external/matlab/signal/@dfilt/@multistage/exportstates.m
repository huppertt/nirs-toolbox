function zf = exportstates(Hd)
%EXPORTSTATES Export the final conditions.

% This should be a private method

%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

zf = {};
for i=1:nstages(Hd),
    zf{i} = exportstates(Hd.Stage(i));
end
