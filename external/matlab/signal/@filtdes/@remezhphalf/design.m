function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Set up design params
N = get(d,'order');

% Get passband frequency, it has been prenormalized
Fpass = get(d,'Fpass');

if Fpass <= 0.5,
    error(message('signal:filtdes:remezhphalf:design:InvalidRange'));
end

b = firhalfband(N,1-Fpass,'high');

% Construct object
Hd = dfilt.dffir(b);



