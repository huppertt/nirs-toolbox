function h = lpcheby2
%LPCHEBY2  Constructor for the lowpass chebyshev type II filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.lpcheby2;

% Call the super's constructor
filterType_construct(h);





