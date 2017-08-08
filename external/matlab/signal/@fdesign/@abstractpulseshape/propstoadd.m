function p = propstoadd(this)
%PROPSTOADD   Return the properties to add to the parent object.

%   Copyright 2008-2009 The MathWorks, Inc.

p = propstoadd(this.CurrentSpecs);

p = [{'Description', 'SamplesPerSymbol'}, p];

% [EOF]
