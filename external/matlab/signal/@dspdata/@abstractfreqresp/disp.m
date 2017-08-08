function disp(this)
%DISP   Display method.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

proplist = reorderprops(this);
snew = reorderstructure(get(this),proplist{:});

val = 'false';
if this.NormalizedFrequency,
    val = 'true';
end
snew.NormalizedFrequency = val;
snew = changedisplay(snew,'NormalizedFrequency',val);
disp(snew);

% [EOF]
