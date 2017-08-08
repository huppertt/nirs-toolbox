function a = getindex(hObj)
%GETINDEX Returns a vector of indexes

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

[dindx, qindx] = getfiltindx(hObj);

a = [qindx; qindx];
a = reshape(a,1,length(qindx)*2);
a = [a dindx];

% [EOF]
