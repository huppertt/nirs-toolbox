function refscalevalues = setrefscalevalues(Hd, refscalevalues)
%SETREFSCALEVALUES   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

validaterefcoeffs(Hd.filterquantizer, 'ScaleValues', refscalevalues);


% [EOF]
