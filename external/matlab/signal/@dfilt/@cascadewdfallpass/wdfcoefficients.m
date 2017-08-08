function coeffs = wdfcoefficients(this)
%WDFCOEFFICIENTS   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

coeffs = this.AllpassCoefficients;
fn = fieldnames(coeffs);
for k = 1:length(fn),
    coeffs = setfield(coeffs,fn{k},convertcoeffs(getfield(coeffs,fn{k})));
end

%--------------------------------------------------------------------------
function gamma = convertcoeffs(c)

switch length(c),
    case 0,
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
