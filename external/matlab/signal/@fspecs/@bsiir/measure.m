function measurements = measure(this, hd, hm)
%MEASURE   Get the measurements.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

measurements = get_bs_measurements(this, hd, hm, [], [], ...
    this.Fpass1, this.Fstop1, this.Fstop2, this.Fpass2, [], [], []);

% [EOF]


