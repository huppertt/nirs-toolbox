function reflattice = setreflattice(Hd, reflattice)
%SETREFLATTICE   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

validaterefcoeffs(Hd.filterquantizer, 'Lattice', reflattice);

% [EOF]
