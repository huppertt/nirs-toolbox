function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

args = getarguments(h,d);

b = firls(get(d, 'Order'), args{:});

% Construct object
Hd = dfilt.dffir(b);

% [EOF]
