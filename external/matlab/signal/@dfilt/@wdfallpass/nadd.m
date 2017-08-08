function N = nadd(this)
%NADD   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

cstruct = wdfcoefficients(reffilter(this));

c = cstruct.Section1;

N = 0;
for k = 1:length(c),
    switch(c(k)),
        case {1,-1},
            N = N + 2;
        case 0,
            % Do nothing
        otherwise
            N = N + 3;
    end
end

% [EOF]
