function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

args = getarguments(h, d);

N = get(d, 'Order');

if rem(N,2),
    % Cannot design odd order bandstops
    error(message('signal:filtdes:firlsbspassstop:design:MustBeEven'));
end

b = firls(N, args{:});

% Construct object
Hd = dfilt.dffir(b);

% [EOF]
