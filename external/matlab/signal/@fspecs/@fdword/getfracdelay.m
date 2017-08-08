function fracdelay = getfracdelay(this)
%GETFRACDELAY   Get the fracdelay.

%   Author(s): V. Pellissier
%   Copyright 2006-2011 The MathWorks, Inc.
%     

fracdelay = this.privFracDelay;
if fracdelay<0,
    error(message('signal:fspecs:fdword:getfracdelay:NonPositiveFracDelay'));
end

if this.NormalizedFrequency,
    Fs = 1;
else
    Fs = this.Fs;
end
if fracdelay>=1/Fs,
    error(message('signal:fspecs:fdword:getfracdelay:InvalidFracDelay', num2str( 1/Fs )));
end

% [EOF]
