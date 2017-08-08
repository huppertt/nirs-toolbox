function varargout = thissfcnparams(Hd)
%THISSFCNPARAMS Returns the parameters for SDSPFILTER

% Author(s): J. Schickler
% Copyright 1988-2002 The MathWorks, Inc.

% Convert the Numerator to a string
num = sprintf('%.25g, ', Hd.Numerator);
num(end-1:end) = [];

varargout = {1, dfobjsfcnparams(Hd), num, '', '0'};

% [EOF]
