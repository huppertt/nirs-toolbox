function reset(Hd)
%RESET Reset the filter.


%   Author: P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

Hd.NumSamplesProcessed = 0;

for k=1:length(Hd.Stage)
    reset(Hd.Stage(k));
end
