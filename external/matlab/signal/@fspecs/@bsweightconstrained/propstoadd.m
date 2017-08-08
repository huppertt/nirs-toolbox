function p = propstoadd(this)
%PROPSTOADD Return the properties to add to the parent object.

%   Copyright 2011 The MathWorks, Inc.

p = {'NormalizedFrequency', 'Fs','FilterOrder','Fpass1',...
  'Fstop1','Fstop2','Fpass2','Passband1Constrained',...
  'StopbandConstrained','Passband2Constrained'};

if this.Passband1Constrained
  p{end+1} = 'Apass1';
end
if this.StopbandConstrained
  p{end+1} = 'Astop';
end
if this.Passband2Constrained
  p{end+1} = 'Apass2';
end

% [EOF]
