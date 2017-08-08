function p = thisprops2add(~,varargin)
%THISPROPS2ADD Return the properties to add to the parent object.

%   Copyright 2011 The MathWorks, Inc.

if ischar(varargin{end})
  varargin(end) = [];
end

idx = strcmp('Fs',varargin);
varargin(idx) = [];

p = {'NormalizedFrequency', 'Fs','FilterOrder','Fstop1',...
  'Fpass1','Fpass2','Fstop2'};

logicalIndex = [];
if length(varargin) > 5
  for idx = 6:length(varargin)
    if islogical(varargin{idx})
      logicalIndex = [logicalIndex idx]; %#ok<AGROW>
    end
  end
end

if ~isempty(logicalIndex)
  Stopband1Constrained = varargin{logicalIndex(1)};
  p{end+1} = 'Stopband1Constrained';
  if Stopband1Constrained && length(varargin)>logicalIndex(1) && ~islogical(varargin{logicalIndex(1)+1})
    p{end+1} = 'Astop1';    
  end
  if length(logicalIndex) > 1
    PassbandConstrained = varargin{logicalIndex(2)};
     p{end+1} = 'PassbandConstrained';
    if PassbandConstrained && length(varargin)>logicalIndex(2) && ~islogical(varargin{logicalIndex(2)+1})
      p{end+1} = 'Apass';    
    end
    if length(logicalIndex) > 2
      Stopband2Constrained = varargin{logicalIndex(3)};
       p{end+1} = 'Stopband2Constrained';
      if Stopband2Constrained && length(varargin)>logicalIndex(3) && ~islogical(varargin{logicalIndex(3)+1})
        p{end+1} = 'Astop2';    
      end
    end
  end
end
    
% Remove the NormalizedFrequency and Fs properties.
p(1:2) = [];

% [EOF]
