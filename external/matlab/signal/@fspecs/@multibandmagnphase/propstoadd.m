function p = propstoadd(this)
%PROPSTOADD   Return the properties to add to the parent object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

p = [{'NormalizedFrequency', 'Fs'}, orderprop(this), {'NBands'}];
for i=1:this.NBands,
    p = [p {sprintf('%s%d%s','B',i,'Frequencies'), ...
        sprintf('%s%d%s','B',i,'FreqResponse')}];
end

% [EOF]
