function h = bsbutter
%BSBUTTER  Constructor for the bandstop butterworth filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.bsbutter;

% Call the super's constructor
filterType_construct(h);

