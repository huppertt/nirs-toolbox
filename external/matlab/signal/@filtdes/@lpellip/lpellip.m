function h = lpellip
%LPELLIP  Constructor for the lowpass elliptic filter type.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.lpellip;

% Call the super's constructor
filterType_construct(h);





