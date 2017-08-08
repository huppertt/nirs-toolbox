function refcoeffs = setrefcoeffs(this, refcoeffs)
%SETREFCOEFFS   Set the refcoeffs.

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

validaterefcoeffs(this.filterquantizer, 'Coefficients', refcoeffs);

% [EOF]
