function specs = thisgetspecs(this)
%THISGETSPECS Get the specs.

%   Copyright 2011 The MathWorks, Inc.

specs.Fpass           = this.Fcutoff;
specs.Fstop           = this.Fcutoff;
specs.Apass           = this.Apass;
specs.Astop           = this.Astop;
specs.FrequencyFactor = this.FrequencyFactor;
specs.Power           = this.Power;

% [EOF]
