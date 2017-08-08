function p = propstoadd(this)
%PROPSTOADD   Return the properties to add to the parent object.

%   Copyright 2008 The MathWorks, Inc.

%The first two properties are added here since they will be removed by the
%inherited thisprops2add methd
p = {'NormalizedFrequency','Fs','NumPeaksOrNotches',...
    'BW','GBW','ShelvingFilterOrder'};

% [EOF]
