function p = propstoadd(this)
%PROPSTOADD   

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

p = fieldnames(this);

% Remove the ResponseType
p(1) = [];

% Remove privFracdelay
p(end) = [];



% [EOF]
