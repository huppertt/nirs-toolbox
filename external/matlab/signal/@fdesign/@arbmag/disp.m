function disp(this)
%DISP Display the design object.

%   Copyright 2011 The MathWorks, Inc.

s = get(this);
s = reorderstructure(s, 'Response', 'Specification', 'Description');

if isfield(s,'Ripple')
  ripple = this.Ripple;
  s = rmfield(s,'Ripple');
  s.Ripple = ripple;
elseif isfield(s,'B1Constrained')
  propNames = {'Response', 'Specification', 'Description',...
    'NormalizedFrequency','Fs','FilterOrder','NBands'};
  for idx = 1:this.NBands
    newCell = {sprintf('%s%d%s','B',idx,'Frequencies')...
      sprintf('%s%d%s','B',idx,'Amplitudes'),...
      sprintf('%s%d%s','B',idx,'Constrained')};
    if this.(sprintf('%s%d%s','B',idx,'Constrained'))
      newCell = [newCell {sprintf('%s%d%s','B',idx,'Ripple')}]; %#ok<*AGROW>
    end
    propNames = [propNames newCell];
  end
  s = reorderstructure(s, propNames{:});
end

if s.NormalizedFrequency
  s = rmfield(s, 'Fs');
end
siguddutils('dispstr', s);

% [EOF]
