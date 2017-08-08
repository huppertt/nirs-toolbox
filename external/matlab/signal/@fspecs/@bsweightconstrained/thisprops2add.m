function p = thisprops2add(~,varargin)
%THISPROPS2ADD Return the properties to add to the parent object.

%   Copyright 2011 The MathWorks, Inc.

if ischar(varargin{end})
  varargin(end) = [];
end

idx = strcmp('Fs',varargin);
varargin(idx) = [];

p = {'NormalizedFrequency', 'Fs','FilterOrder','Fpass1',...
  'Fstop1','Fstop2','Fpass2'};

logicalIndex = [];
if length(varargin) > 5
  for idx = 6:length(varargin)
    if islogical(varargin{idx})
      logicalIndex = [logicalIndex idx]; %#ok<AGROW>
    end
  end
end

if ~isempty(logicalIndex)
  Passband1Constrained = varargin{logicalIndex(1)};
  p{end+1} = 'Passband1Constrained';
  if Passband1Constrained && length(varargin)>logicalIndex(1) && ~islogical(varargin{logicalIndex(1)+1})
    p{end+1} = 'Apass1';    
  end
  if length(logicalIndex) > 1
    StopbandConstrained = varargin{logicalIndex(2)};
     p{end+1} = 'StopbandConstrained';
    if StopbandConstrained && length(varargin)>logicalIndex(2) && ~islogical(varargin{logicalIndex(2)+1})
      p{end+1} = 'Astop';    
    end
    if length(logicalIndex) > 2
      Passband2Constrained = varargin{logicalIndex(3)};
       p{end+1} = 'Passband2Constrained';
      if Passband2Constrained && length(varargin)>logicalIndex(3) && ~islogical(varargin{logicalIndex(3)+1})
        p{end+1} = 'Apass2';    
      end
    end
  end
end
    
% Remove the NormalizedFrequency and Fs properties.
p(1:2) = [];

% [EOF]
