function Hd = design(h,d)
%Design  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Set up design params
N = get(d,'order');
Fc1 = get(d,'Fc1');
Fc2 = get(d,'Fc2');
Fc = [Fc1, Fc2];

win = generatewindow(d);

scaleflag = determinescaleflag(d);


b = fir1(N,Fc,'bandpass',win,scaleflag);

% Construct object
Hd = dfilt.dffir(b);

