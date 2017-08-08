function g = thisnormalize(Hd)
%THISNORMALIZE   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

lat = Hd.reflattice;
g = max(abs(lat));
Hd.reflattice= lat/g;

% [EOF]
