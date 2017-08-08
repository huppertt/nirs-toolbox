function bands = set_bands(this, bands)
%SET_BANDS   PreSet function for the 'bands' property.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

if ~isnumeric(bands) || round(bands)~=bands || bands<1 || bands>10,
    error(message('signal:fspecs:abstractmultiband:set_bands:InvalidNumberofBands'));
end

% [EOF]
