function p = propstoadd(this)
%PROPSTOADD   Return the properties to add to the parent object.

%   Copyright 2008 The MathWorks, Inc.

p = fieldnames(this);
p = {p{1:4},p{6},p{8},p{7},p{5}};

p(1) = []; % All but the responsetype.

% [EOF]
