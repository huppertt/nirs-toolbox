function h = hpcheby1
%HPCHEBY1  Constructor for the highpass chebyshev type I filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.hpcheby1;

% Call the super's constructor
filterType_construct(h);




