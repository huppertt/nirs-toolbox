function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

L = get(d, 'InterpolationFactor');

[Fpass, Fstop, Dpass, Dstop] = getdesignspecs(h,d);

[h, g, b] = ifir(L, 'low', [Fpass Fstop], [Dpass Dstop], get(d, 'Optimization'));

% Construct object
Hd = cascade(dfilt.dffir(h), dfilt.dffir(g));

if L*(1-Fpass) < 1,
    Hd = parallel(Hd, dfilt.dffir(b));
end

% [EOF]
