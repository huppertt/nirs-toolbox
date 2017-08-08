function varargout = thissfcnparams(Hd)
%THISSFCNPARAMS Returns the parameters for SDSPFILTER

% Author(s): J. Schickler
% Copyright 1988-2002 The MathWorks, Inc.

num = Hd.Numerator;
den = Hd.Denominator;

[filtertype, filterstrt, ic] = dfobjsfcnparams(Hd);

% Normalize the coefficients
if den(1) ~= 1,
    num = num ./ den(1);
    den = den ./ den(1);
else
    
    % Check for the all-pole case
    if length(num) == 1 & num == 1,
        filtertype = 2;
        filterstrt = 1;
    end
end

% Convert the Numerator and Denominator to strings
num = sprintf('%.25g, ', num);
num(end-1:end) = [];
den = sprintf('%.25g, ', den);
den(end-1:end) = [];

varargout = {filtertype, filterstrt, num, den, ic};

% [EOF]
