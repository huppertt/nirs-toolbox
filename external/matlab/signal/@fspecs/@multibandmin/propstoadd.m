function p = propstoadd(this)
%PROPSTOADD Return the properties to add to the parent object.

%   Copyright 2011 The MathWorks, Inc.

p = [{'NormalizedFrequency', 'Fs'}, {'NBands'}];
for i=1:this.NBands,
    p = [p {sprintf('%s%d%s','B',i,'Frequencies'), ...
        sprintf('%s%d%s','B',i,'Amplitudes'), ...
        sprintf('%s%d%s','B',i,'Ripple')}]; %#ok<*AGROW>
end

% [EOF]
