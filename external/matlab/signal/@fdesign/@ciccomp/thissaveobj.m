function s = thissaveobj(this)
%THISSAVEOBJ Save this object.

%   Copyright 2005-2011 The MathWorks, Inc.

s.DifferentialDelay   = this.DifferentialDelay;
s.NumberOfSections    = this.NumberOfSections;
s.CICRateChangeFactor = this.CICRateChangeFactor; % Property added in R2011b

% [EOF]