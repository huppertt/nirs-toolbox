function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Set up design params
N = get(d,'order');

[F,A,W] = getNumericSpecs(h,d);

b = firls(N,F,A,W,'hilbert');

% Construct object
Hd = dfilt.dffir(b);



