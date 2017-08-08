function refladder = setrefladder(Hd, refladder)
%SETREFLADDER   

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

validaterefcoeffs(Hd.filterquantizer, 'Ladder', refladder);

% [EOF]
