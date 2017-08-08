function specs = thisgetspecs(this)
%THISGETSPECS Get the specs.

%   Copyright 2005-2011 The MathWorks, Inc.

specs.Fpass               = this.Fstop;
specs.Fstop               = this.Fstop;
specs.Apass               = this.Apass;
specs.Astop               = this.Astop;
specs.FrequencyFactor     = this.FrequencyFactor;
specs.Power               = this.Power;
specs.CICRateChangeFactor = this.CICRateChangeFactor;


% [EOF]
