function h = lpnorm
%LPNORM  Constructor for the lpnorm specifications object.
%
%   Outputs:
%       h - Handle to the lpnorm object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


h = filtdes.lpnorm;
% Call the real constructor
lpnorm_construct(h);
