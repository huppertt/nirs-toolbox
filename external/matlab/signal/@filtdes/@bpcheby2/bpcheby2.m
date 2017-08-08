function h = bpcheby2
%BPCHEBY2  Constructor for the bandpass chebyshev type II filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.bpcheby2;

% Call the super's constructor
filterType_construct(h);





