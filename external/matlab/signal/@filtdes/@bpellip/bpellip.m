function h = bpellip
%BPELLIP  Constructor for the bandpass elliptic filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.bpellip;

% Call the super's constructor
filterType_construct(h);





