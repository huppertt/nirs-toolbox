function p = propstoadd(this)
%PROPSTOADD   

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

p = fieldnames(this);
p(1) = []; % All but the responsetype.

% [EOF]
