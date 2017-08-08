function denormalize(Hd)
%DENORMALIZE Undo normalization applied by NORMALIZE.
% 
%   See also: NORMALIZE

%   Author(s): V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

for k=1:length(Hd.Stage)
    denormalize(Hd.Stage(k));
end


% [EOF]
