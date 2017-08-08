function loadreferencecoefficients(this, s)
%LOADREFERENCECOEFFICIENTS   Load the reference coefficients.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if s.version.number < 2
    set(this, 'Lattice', s.Lattice, ...
        'Ladder', s.Ladder);
else
    set(this, 'Lattice', s.reflattice, ...
        'Ladder', s.refladder);
end

% [EOF]
