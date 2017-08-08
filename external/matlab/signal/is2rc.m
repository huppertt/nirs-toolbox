function k = is2rc(inv_sin)
%IS2RC  Convert inverse sine parameters to reflection coefficients.
%   K = IS2RC(INV_SIN) returns the reflection coefficients corresponding 
%   to the inverse sine parameters, INV_SIN. 
%
%   % Example:
%   %   Convert the following inverse sine parameters to reflection 
%   %   coefficients.
%   %   isin = [0.2000 0.8727 0.0020 0.0052 -0.0052];
%
%   isin = [0.2000 0.8727 0.0020 0.0052 -0.0052];
%   k = is2rc(isin)     % Gives Reflection coefficients
%
%   See also RC2IS, POLY2RC, AC2RC, LAR2RC.

%   Reference: J.R. Deller, J.G. Proakis, J.H.L. Hansen, "Discrete-Time 
%   Processing of Speech Signals", Prentice Hall, Section 7.4.5.
%
%   Author(s): A. Ramasubramanian
%   Copyright 1988-2002 The MathWorks, Inc.

if ~isreal(inv_sin),
    error(message('signal:is2rc:MustBeReal'));
end

k = sin(inv_sin*pi/2);

% [EOF] is2rc.m
