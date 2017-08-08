function coeffs = set_coeffs(this, coeffs)
%SET_COEFFS   PreSet function for the 'coeffs' property.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

% Convert struct to cell
coeffs = struct2cell(coeffs);

% Check that coeffs are valid
validate_coeffs(this,coeffs);

% Make sure to clear metadata
clearmetadata(this);


oldnsecs = length(this.refallpasscoeffs);
oldncoeffs = zeros(oldnsecs,1);
for k = 1:oldnsecs,
    oldncoeffs(k) = length(this.refallpasscoeffs{k});
end

% Set the reference coefficients
this.refallpasscoeffs = coeffs;

% Quantize the coefficients
quantizecoeffs(this);

% If number of coeffs changes, flush states
nsecs = length(coeffs);
if  oldnsecs ~= nsecs,
    % Number of sections has changed
    reset(this);
else
    ncoeffs = zeros(nsecs,1);
    for k = 1:nsecs,
        ncoeffs(k) = length(coeffs{k});
    end
    if any(oldncoeffs ~= ncoeffs),
        % Length of a section has changed
        reset(this);
    end
end

% Hold an empty to not duplicate storage
coeffs = [];

% [EOF]
