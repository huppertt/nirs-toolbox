function thisunnormalize(Hd, g)
%THISUNNORMALIZE   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

lat = Hd.reflattice;
Hd.reflattice= lat*g;

% [EOF]
