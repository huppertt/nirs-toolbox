function thiscopy(this, hOldObject)
%THISCOPY Copy properties specific to the fdesign.ciccomp class.

%   Copyright 2005-2011 The MathWorks, Inc.

set(this, 'DifferentialDelay', hOldObject.DifferentialDelay, ...
    'NumberOfSections', hOldObject.NumberOfSections,...
    'CICRateChangeFactor', hOldObject.CICRateChangeFactor );

% [EOF]
