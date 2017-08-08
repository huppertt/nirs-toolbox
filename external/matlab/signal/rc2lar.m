function g = rc2lar(k)
%RC2LAR Convert reflection coefficients to log area ratios.
%   G = RC2LAR(K) returns the log area ratios corresponding to the reflection
%   coefficients, K.
%
%   The log area ratio is defined by G = log((1+k)/(1-k)) , where the K is the
%   reflection coefficient.
%
%   % Example:
%   %   Convert the following reflection coefficients into log area ratio
%   %   parameters.
%   %   k = [0.3090    0.9800    0.0031    0.0082   -0.0082]
%
%   k = [0.3090 0.9801 0.0031 0.0082 -0.0082];
%   g = rc2lar(k)       % Gives log area ratios
%
%   See also LAR2RC, RC2POLY, RC2AC, RC2IS.

%   References:
%   [1] J. Makhoul, "Linear Prediction: A Tutorial Review," Proc. IEEE,
%   Vol.63, No.4, pp.561-580, Apr 1975.
%   [2] ITU-T Recommendation G.729 (03/96)
%
%   Author(s): A. Ramasubramanian
%   Copyright 1988-2012 The MathWorks, Inc.

if ~isreal(k),
    error(message('signal:rc2lar:NotSupported'));
end

if max(abs(k)) >= 1,
    error(message('signal:rc2lar:InvalidRange'));
end

% Use the relation, atanh(x) = (1/2)*log((1+k)/(1-k))

g = -2*atanh(-k);

% [EOF] rc2lar.m
