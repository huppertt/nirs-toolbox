function thisunnormalize(Hd, g)
%THISUNNORMALIZE   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

B = Hd.refB;
C = Hd.refC;
D = Hd.refD;
Hd.refB = B*g(1);
Hd.refC = C*g(2);
Hd.refD = D*g(3);


% [EOF]
