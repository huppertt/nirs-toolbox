function F0 = set_F0(this, F0)
%SET_F0 PreSet function for the 'F0' property

%   Copyright 2008 The MathWorks, Inc.

if isequal(F0,1) || isequal(F0,0)
    error(message('signal:fspecs:parameqaudioqa:set_F0:invalidSpecs'));
end
% [EOF]
