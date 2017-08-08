function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

L = get(d, 'InterpolationFactor');

[Fstop, Fpass, Dstop, Dpass] = getdesignspecs(h,d);

[h, g, b] = ifir(L, 'high', [Fstop Fpass], [Dstop Dpass], get(d, 'Optimization'));

% Construct object
Hd = cascade(dfilt.dffir(h), dfilt.dffir(g));

if L*Fpass < 1,
    Hd = parallel(Hd, dfilt.dffir(b));
end

% [EOF]
