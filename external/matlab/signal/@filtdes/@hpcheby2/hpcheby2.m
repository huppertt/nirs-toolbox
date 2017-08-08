function h = hpcheby2
%HPCHEBY2  Constructor for the highpass chebyshev type II filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.hpcheby2;

% Call the super's constructor
filterType_construct(h);





