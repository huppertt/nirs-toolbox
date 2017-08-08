function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Set up design params
N = get(d,'order');

win = generatewindow(d);

L = get(d,'band');

b = firnyquist(N,L,win);

% Construct object
Hd = dfilt.dffir(b);



