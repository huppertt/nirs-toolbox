function h = bsellip
%BSELLIP  Constructor for the bandstop elliptic filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.bsellip;

% Call the super's constructor
filterType_construct(h);




