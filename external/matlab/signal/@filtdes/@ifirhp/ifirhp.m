function h = ifirhp
%IFIRHP Constructor for this object.
%
%   Inputs:
%
%   Outputs:
%       h - Handle to this object

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

h = filtdes.ifirhp;

% Call the super's constructor
filterType_construct(h);
