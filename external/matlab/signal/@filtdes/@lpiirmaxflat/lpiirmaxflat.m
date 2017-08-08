function h = lpiirmaxflat
%LPIIRMAXFLAT  Constructor for the lowpass maxflat filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.lpiirmaxflat;

% Call the super's constructor
filterType_construct(h);

