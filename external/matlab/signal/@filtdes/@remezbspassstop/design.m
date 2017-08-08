function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

args = getarguments(h, d);

N = get(d, 'Order');

% Check for valid order
if rem(N, 2),
    error(message('signal:filtdes:remezbspassstop:design:MustBeEven'));
end

dens = get(d,'DensityFactor');

b = remez(N, args{:}, {dens});

% Construct object
Hd = dfilt.dffir(b);

% [EOF]
