function refsosmatrix = setrefsosmatrix(Hd, refsosmatrix)
%SETREFSOSMATRIX   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

validaterefcoeffs(Hd.filterquantizer, 'sosMatrix', refsosmatrix);

% [EOF]
