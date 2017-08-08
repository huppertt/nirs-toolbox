function h = hpbutter
%HPBUTTER  Constructor for the highpass butterworth filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.hpbutter;

% Call the super's constructor
filterType_construct(h);





