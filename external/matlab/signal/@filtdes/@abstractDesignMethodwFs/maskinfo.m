function cmd = maskinfo(d)
%MASKINFO Returns the a cell of structures containing the mask information.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd = maskinfo(get(d, 'responseTypeSpecs'), d);

% [EOF]
