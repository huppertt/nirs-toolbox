function n = nmult(this,optimones,optimnegones)
%NMULT Returns the number of multipliers  

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

% Coefficients
c = coefficients(reffilter(this));
% Structure dependent factor e.g. dfsymfir uses only half the
% coefficients
[f offset] = multfactor(this);

n = 0;

% Number of non-trivial coefficients
for i=1:length(c)
    coeffs = c{i};
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
    n = n + ceil(length(coeffs)*f(i))-offset(i);
end

n = max(n,0);


% [EOF]
