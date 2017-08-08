function specs = thisgetspecs(this)
%THISGETSPECS   

%   Copyright 2008 The MathWorks, Inc.

specs.CombType = this.CombType;
specs.NormalizedFrequency = this.NormalizedFrequency;
specs.Fs = this.Fs;
specs.FilterOrder = this.NumPeaksOrNotches;
specs.BW = this.BW;
specs.PeakNotchFrequencies = this.PeakNotchFrequencies;
% [EOF]
