function p = propstoadd(this)
%PROPSTOADD   Return the properties to add to the parent object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

p = fieldnames(this);

p(1) = []; % All but the responsetype.

% [EOF]
