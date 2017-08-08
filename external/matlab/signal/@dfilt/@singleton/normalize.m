function varargout = normalize(Hd)
%NORMALIZE Normalize coefficients between -1 and 1.
%   G = NORMALIZE(Hd) normalizes the coefficients between -1 and 1 and
%   returns the gain G due to normalization.  Subsequent calls to NORMALIZE
%   will not change the coefficients and G will always return the gain used
%   in the first normalization.
%
%   See also DFILT/DENORMALIZE.

%   Author(s): V. Pellissier
%   Copyright 1988-2005 The MathWorks, Inc.

sv = Hd.privnormGain;
if isempty(sv)
    sv = thisnormalize(Hd);
    Hd.privnormGain = sv;

    % Fire listener to re-quantize coeffs
    quantizecoeffs(Hd);
end

if nargout==1,
    varargout = {sv};
end

% [EOF]
