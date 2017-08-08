function p = propstoadd(this)
%PROPSTOADD Return the properties to add to the parent object.

%   Copyright 2010 The MathWorks, Inc.

p = [{'NormalizedFrequency', 'Fs'}, orderprop(this), {'NBands'}];
for i=1:this.NBands,
    p = [p {sprintf('%s%d%s','B',i,'Frequencies'), ...
        sprintf('%s%d%s','B',i,'GroupDelay')}]; %#ok<AGROW>
end

% [EOF]
