function n = nmult(this,optimones,optimnegones)
%NMULT Returns the number of multipliers  

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

% Coefficients
c = coefficients(this);
coeffs = c{1};
nphases = size(coeffs,2);
coeffs = coeffs(:);
% Exclude zeros
coeffs(find(coeffs==0)) = [];
if optimones,
    % Exclude ones
    coeffs(find(coeffs==1)) = [];
end
if optimnegones,
    % Exclude neg ones
    coeffs(find(coeffs==-1)) = [];
end
% Apply structure dependent factor e.g. dfsymfir uses only half the
% coefficients
n = nphases-1 + length(coeffs);

n = max(n,0);


% [EOF]
