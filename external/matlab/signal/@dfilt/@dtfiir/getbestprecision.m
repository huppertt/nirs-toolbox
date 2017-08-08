function s = getbestprecision(h)
%GETBESTPRECISION Return best precision for Product and Accumulator

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

s = getbestprecision(h.filterquantizer);

% [EOF]
