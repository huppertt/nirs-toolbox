function denormalize(Hd)
%DENORMALIZE   Reverse coefficient changes made by NORMALIZE.
%   DENORMALIZE(Hd) reverses coefficient changes made by NORMALIZE.  This
%   method will not change the coefficients if called before NORMALIZE or
%   on subsequent calls.
% 
%   See also DFILT/NORMALIZE.

%   Author(s): V. Pellissier
%   Copyright 1988-2005 The MathWorks, Inc.

sv = Hd.privnormGain;
if isempty(sv),
    % No op
    return; 
end

thisunnormalize(Hd,sv);
Hd.privnormGain = [];
% Fire listener to re-quantize coeffs
quantizecoeffs(Hd);

% [EOF]
