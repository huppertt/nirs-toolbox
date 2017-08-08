function Hd = design(h,d)
%DESIGN  Method to design the filter given the specs.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

[Fpass, Dpass] = getdesignspecs(h, d);

if Fpass <= 0.5,
    error(message('signal:filtdes:remezhphalfmin:design:InvalidRange'));
end

b = firhalfband('minorder',1-Fpass,Dpass,'high');

% Construct object
Hd = dfilt.dffir(b);



