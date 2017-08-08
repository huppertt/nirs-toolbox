function d = iirlpnormc
%IIRLPNORMC  Constructor for this design method object.
%
%   Outputs:
%       d - Handle to the design method object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.



d = filtdes.iirlpnormc;

% Call the super constructor
iirlpnorm_construct(d);

% Overwrite the tag
set(d,'Tag','IIR Constrained least P-th norm');








