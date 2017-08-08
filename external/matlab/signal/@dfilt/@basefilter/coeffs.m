function s = coeffs(this)
%COEFFS   Return the coefficients in a structure.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

N = length(this);
fn = coefficientnames(this);

if N==1,
    for j = 1:length(fn),
        s.(fn{j}) = get(this, fn{j});
    end
else
    % Build a structure with field names equal to the coefficient names.
    for indx = 1:N
        aux = fn{indx};
        for j = 1:length(aux),
            ss.(aux{j}) = get(this(indx), aux{j});
        end
        s{indx} = ss;
        clear ss;
    end
end

% [EOF]
