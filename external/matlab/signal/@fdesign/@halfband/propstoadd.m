function p = propstoadd(this)
%PROPSTOADD   Return the properties to add to the parent object.

%   Copyright 2007 The MathWorks, Inc.

p = propstoadd(this.CurrentSpecs);

p = {'Description','Type', p{:}};

% [EOF]
