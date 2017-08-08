function p = propstoadd(this)
%PROPSTOADD Return the properties to add to the parent object.

%   Copyright 2011 The MathWorks, Inc.

p = {'NormalizedFrequency', 'Fs','FilterOrder','NBands'};
for i=1:this.NBands,
  p = [p {sprintf('%s%d%s','B',i,'Frequencies'), ...
    sprintf('%s%d%s','B',i,'Amplitudes'), ...
    sprintf('%s%d%s','B',i,'Constrained')}]; %#ok<*AGROW>
  
  % Do not add ripple properties if band has not been specified as
  % constrained. 
  if this.(sprintf('%s%d%s','B',i,'Constrained'))
    p = [p {sprintf('%s%d%s','B',i,'Ripple')}];
  end
end

% [EOF]
