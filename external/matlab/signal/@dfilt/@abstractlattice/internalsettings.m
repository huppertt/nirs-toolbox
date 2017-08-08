function s = internalsettings(h)
%INTERNALSETTINGS Returns the fixed-point settings viewed by the algorithm.  

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

s = internalsettings(h.filterquantizer);

% [EOF]
