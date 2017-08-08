function wt = set_weightingtype(this, wt) 
%SET_WEIGHTINGTYPE PreSet function for the 'WeightingType' property

%   Copyright 2009 The MathWorks, Inc.

if strcmpi(wt,'itur4684');
    %The default 48e3 sampling frequency is not enough to cover the entire range
    %of 64 KHz specified by the ITU-R 468-4 standard so we set the default Fs in
    %this case to 80 KHz (64 KHz is not a typical audio sampling frequency but
    %80 KHz is). We do this to ensure that the default design meets the specs. 
    this.DefaultFs = 80e3;
else
    this.DefaultFs = 48e3;
end
    
% [EOF]
