function coeffs = get_coeffs(this, coeffs)
%GET_COEFFS   PreGet function for the 'coeffs' property.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

% Use getnumerator since it simply casts to double
for k = 1:length(this.privallpasscoeffs),
    coeffs{k} = getnumerator(this.filterquantizer, this.privallpasscoeffs{k});
    flds{k}      = sprintf('Section%d',k);
end
if ~isempty(coeffs),
    coeffs = cell2struct(coeffs,flds,2);
end

% [EOF]
