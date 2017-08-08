function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

specs.Frequencies = this.Frequencies;
specs.Magnitudes = abs(this.FreqResponse);
specs.Phases = angle(this.FreqResponse);


% [EOF]
