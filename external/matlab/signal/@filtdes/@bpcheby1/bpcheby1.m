function h = bpcheby1
%BPCHEBY1  Constructor for the bandpass chebyshev type I filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.bpcheby1;

% Call the super's constructor
filterType_construct(h);




