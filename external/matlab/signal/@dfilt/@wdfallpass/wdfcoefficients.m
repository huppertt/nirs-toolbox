function c = wdfcoefficients(this,coeffs)
%WDFCOEFFICIENTS   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

if nargin < 2,
    c.Section1 = convertcoeffs(this.AllpassCoefficients);
else
    c.Section1 = convertcoeffs(coeffs);
end


%--------------------------------------------------------------------------
function gamma = convertcoeffs(c)

switch length(c),
    case 0, 
        % Empty coeffs
        gamma = [];
    case 1,
        gamma = c;
    case 2,
        gamma = [c(2) c(1)/(1+c(2))];
    case 4,
        % This works only when c is of the form [0 c(2) 0 c(4)], i.e.
        % halfband linear phase
        gamma = [c(4) 0 c(2)/(1+c(4)) 0];
end

% [EOF]
