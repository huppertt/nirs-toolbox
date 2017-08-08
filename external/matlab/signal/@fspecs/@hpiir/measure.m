function measurements = measure(this, hd, hm)
%MEASURE   Get the measurements.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

measurements = get_hp_measurements(this, hd, hm, [], this.Fstop, this.Fpass, [], []);

% [EOF]
