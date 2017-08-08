function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

[Fpass, Dpass] = getdesignspecs(h, d);

b = firhalfband('minorder',1-Fpass,Dpass,'kaiser','high');

% Construct object
Hd = dfilt.dffir(b);

% [EOF]
