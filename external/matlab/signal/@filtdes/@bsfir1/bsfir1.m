function h = bsfir1
%BSFIR1  Constructor for this object.
%
%   Outputs:
%       h - Handle to this object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


h = filtdes.bsfir1;

% Call the super's constructor
filterType_construct(h);

