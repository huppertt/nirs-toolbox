function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Set up design params
N = get(d,'Order');

% Get frequency specs, they have been prenormalized
Fc = get(d,'Fc');

[b,a,b1,b2,s,g] = maxflat(N,'sym',Fc);

% Construct object
Hd = dfilt.df2sos(s,g);



