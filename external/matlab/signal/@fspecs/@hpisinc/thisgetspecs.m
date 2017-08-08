function specs = thisgetspecs(this)
%THISGETSPECS Get the specs.

%   Copyright 2011 The MathWorks, Inc.

specs.Fpass           = this.Fpass;
specs.Fstop           = this.Fstop;
specs.Apass           = NaN;
specs.Astop           = NaN;
specs.FrequencyFactor = this.FrequencyFactor;
specs.Power           = this.Power;

% [EOF]
