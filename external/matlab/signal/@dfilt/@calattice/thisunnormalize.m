function thisunnormalize(Hd, g)
%THISUNNORMALIZE   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

lat1 = Hd.refAllpass1;
lat2 = Hd.refAllpass2;
Hd.refAllpass1 = lat1*g(1);
Hd.refAllpass2 = lat2*g(2);


% [EOF]
