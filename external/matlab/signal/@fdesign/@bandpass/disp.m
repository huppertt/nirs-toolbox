function disp(this)
%DISP Display the design object.

%   Copyright 2011 The MathWorks, Inc.

s = get(this);
s = reorderstructure(s, 'Response', 'Specification', 'Description');

if isfield(s,'Stopband1Constrained')
  propNames = {'Response', 'Specification', 'Description',...
    'NormalizedFrequency','Fs','FilterOrder','Fstop1','Fpass1',...
    'Fpass2','Fstop2'};
  
  propNames{end+1} = 'Stopband1Constrained';
  if this.Stopband1Constrained
    propNames{end+1} = 'Astop1';
  end
  propNames{end+1} = 'PassbandConstrained';
  if this.PassbandConstrained
    propNames{end+1} = 'Apass';
  end  
  propNames{end+1} = 'Stopband2Constrained';
  if this.Stopband2Constrained
    propNames{end+1} = 'Astop2';
  end  
  s = reorderstructure(s, propNames{:});
end

if s.NormalizedFrequency
  s = rmfield(s, 'Fs');
end
siguddutils('dispstr', s);

% [EOF]
