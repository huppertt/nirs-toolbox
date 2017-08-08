function newfq = set_filterquantizer(this, newfq)
%SET_FILTERQUANTIZER   PreSet function for the 'filterquantizer' property.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

oldfq = get(this, 'filterquantizer');

set(this, 'privfilterquantizer', newfq);

% Remove current dynamic properties
rmprops(this, oldfq);

set_ncoeffs(newfq, [naddp1(this) size(this.Coefficients,2)]);

try
    % Quantize the coefficients
    quantizecoeffs(this);
    
    % Quantize the fracdelay
    quantizefd(this);
catch ME
    
    % Get the old quantizer because the new one errors.
    set(this, 'privfilterquantizer', oldfq);

    newfq = oldfq;
    quantizecoeffs(this);
    quantizefd(this);

    addprops(this, newfq);
    throwAsCaller(ME);
end

addprops(this, newfq);

validatestates(this);

% Store nothing, its stored in private.
newfq = [];

% [EOF]
