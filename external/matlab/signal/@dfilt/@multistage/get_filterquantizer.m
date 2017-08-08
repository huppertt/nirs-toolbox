function fq = get_filterquantizer(this, fq)
%get the filterquantizer for the first stage of multistage filters.

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

fq = get_filterquantizer(this.Stage(1));

% [EOF]
