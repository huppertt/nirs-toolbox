function p = propstoadd(this)
%PROPSTOADD   Return the properties to add to the parent object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

p = fieldnames(this);
p = {p{1:3},p{6:7},p{10},p{4:5},p{8:9},p{11:end}};

p(1) = []; % All but the responsetype.

% [EOF]
