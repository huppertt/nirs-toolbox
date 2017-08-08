function fr = whichframes(this)
%WHICHFRAMES   Returns the frames for the classic iir.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

fr = dmom_whichframes(this);

if ~isspecify(this)
    fr(end).constructor = 'siggui.classiciiroptsframe';
    fr(end).setops      = {'MatchExactly', get(this, 'MatchExactly')};
end

% [EOF]
