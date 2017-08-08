function h = freqmagweight
%freqmagweight  Constructor for this object.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


h = filtdes.freqmagweight;

% Call alternate constructor
fmw_construct(h);

