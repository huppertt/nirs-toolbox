function d = firlpnorm
%FIRLPNORM  Constructor for this design method.
%
%   Outputs:
%       d - Handle to the design method object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

d = filtdes.firlpnorm;

% Call super's constructor
lpnorm_construct(d);

% Set the tag
set(d,'Tag','FIR least P-th norm');







