function reset(this)
%RESET   Reset the object.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

p = propstocopy(this);

p = [{'Specification'}, p];

c = get(this, 'CapturedState');

for indx = 1:length(p)
    set(this, p{indx}, c.(p{indx}));
end

allSpecs = get(this, 'AllSpecs');

for indx = 1:length(allSpecs)
    f = strrep(class(allSpecs(indx)), '.', '_');
    setstate(allSpecs(indx), c.(f));
end

% [EOF]
