function h = bscheby2
%BSCHEBY2  Constructor for the bandstop chebyshev type II filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.bscheby2;

% Call the super's constructor
filterType_construct(h);





