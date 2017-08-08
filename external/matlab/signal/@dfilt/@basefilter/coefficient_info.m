function [p, v] = coefficient_info(this)
%COEFFICIENT_INFO   Get the coefficient information for this filter.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

coeffs = coefficients(this);
if length(coeffs) == 1
    p = {getString(message('signal:dfilt:info:FilterLength'))};
    v = {sprintf('%d', length(coeffs{1}))};
else
    coeffnames = coefficientnames(this);
    for indx = 1:length(coeffs)
        p{indx} = sprintf('%s %s', coeffnames{indx}, ...
                          getString(message('signal:dfilt:info:Length')));
        v{indx} = sprintf('%d', length(coeffs{indx}));
    end
end

% [EOF]
