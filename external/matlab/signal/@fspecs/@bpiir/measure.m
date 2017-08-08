function measurements = measure(this, hd, hm)
%MEASURE   Get the measurements.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

measurements = get_bp_measurements(this, hd, hm, [], [], ...
    this.Fstop1, this.Fpass1, this.Fpass2, this.Fstop2, [], [], []);

% [EOF]

