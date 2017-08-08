function h = lpcheby1
%LPCHEBY1  Constructor for the lowpass chebyshev type I filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.lpcheby1;

% Call the super's constructor
filterType_construct(h);





