function inv_sin = rc2is(k)
%RC2IS Convert reflection coefficients to inverse sine parameters.
%   INV_SIN = RC2IS(K) returns the inverse sine parameters corresponding 
%   to the reflection coefficients, K.
%
%   % Example:
%   %   Convert the following reflection coefficients to inverse sine 
%   %   parameters.
%   %   k = [0.3090    0.9800    0.0031    0.0082   -0.0082]
%
%   k = [0.3090 0.9801 0.0031 0.0082 -0.0082];
%   isin = rc2is(k)     % Gives inverse sine parameters
%
%   See also IS2RC, RC2POLY, RC2AC, RC2LAR.

%   Reference: J.R. Deller, J.G. Proakis, J.H.L. Hansen, "Discrete-Time 
%   Processing of Speech Signals", Prentice Hall, Section 7.4.5.
%
%   Author(s): A. Ramasubramanian
%   Copyright 1988-2002 The MathWorks, Inc.

if ~isreal(k),
 error(message('signal:rc2is:NotSupported'));
end           

if max(abs(k)) >= 1,
    error(message('signal:rc2is:InvalidRange'));
end

inv_sin = (2/pi)*asin(k);

% [EOF] rc2is.m
