function str = getinfoheader(Hm)
%GETINFOHEADER

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

% Setup the title string
if isreal(Hm), typestr = 'real';
else,         typestr = 'complex'; end

if isfir(Hm), rtype = 'FIR';
else,        rtype = 'IIR'; end

str = sprintf('Discrete-Time %s Farrow Filter (%s)', rtype, typestr);

% [EOF]
