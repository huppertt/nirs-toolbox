function n = nmult(this,optimones,optimnegones)
%NMULT Returns the number of multipliers 

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

% Coefficients
c = struct2cell(wdfcoefficients(reffilter(this)));

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

end

n = length(coeffs);



% [EOF]
