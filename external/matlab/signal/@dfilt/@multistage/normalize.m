function varargout = normalize(Hd)
%NORMALIZE Normalize coefficients of each section between -1 and 1.
%   G = NORMALIZE(Hd) returns the gains due to normalization.

%   See also: DENORMALIZE

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

for k=1:length(Hd.Stage)
    sv{k} = normalize(Hd.Stage(k));
end

if nargout==1,
    varargout = {sv};
end



% [EOF]
