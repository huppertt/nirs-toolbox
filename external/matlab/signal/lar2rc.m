function k = lar2rc(g)
%LAR2RC Convert log area ratios to reflection coefficients.
%   K = LAR2RC(G) returns the reflection coefficients, K, based on the 
%   log area ratios, G.
%
%   % Example:
%   %   Convert the following log area ratio parameters to reflection 
%   %   coefficients
%   %   g = [0.6389    4.5989    0.0063    0.0163   -0.0163];
%
%   g = [0.6389    4.5989    0.0063    0.0163   -0.0163];
%   k = lar2rc(g)       % Reflection coefficients
%
%   See also RC2LAR, POLY2RC, AC2RC, IS2RC.

%   References:
%   [1] J. Makhoul, "Linear Prediction: A Tutorial Review," Proc. IEEE,
%   Vol.63, No.4, pp.561-580, Apr 1975.
%   [2] ITU-T Recommendation G.729 (03/96)
%
%   Author(s): A. Ramasubramanian
%   Copyright 1988-2002 The MathWorks, Inc.

if ~isreal(g),
    error(message('signal:lar2rc:MustBeReal'));
end

% Use the relation, tanh(x) = (1-exp(2x))/(1+exp(2x))
k = -tanh(-g/2);

% [EOF] lar2rc.m
