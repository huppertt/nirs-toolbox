function disp(this)
%DISP Display the design object.

%   Copyright 2011 The MathWorks, Inc.

s = get(this);
s = reorderstructure(s, 'Response', 'Specification', 'Description');

if isfield(s,'Passband1Constrained')
  propNames = {'Response', 'Specification', 'Description',...
    'NormalizedFrequency','Fs','FilterOrder','Fpass1','Fstop1',...
    'Fstop2','Fpass2'};
  
  propNames{end+1} = 'Passband1Constrained';
  if this.Passband1Constrained
    propNames{end+1} = 'Apass1';
  end
  propNames{end+1} = 'StopbandConstrained';
  if this.StopbandConstrained
    propNames{end+1} = 'Astop';
  end  
  propNames{end+1} = 'Passband2Constrained';
  if this.Passband2Constrained
    propNames{end+1} = 'Apass2';
  end  
  s = reorderstructure(s, propNames{:});
end

if s.NormalizedFrequency
  s = rmfield(s, 'Fs');
end
siguddutils('dispstr', s);

% [EOF]
