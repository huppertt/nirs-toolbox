function disp(this)
%DISP   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

s = get(this);
s = reorderstructure(this,s);

if s.NormalizedFrequency,
    nfval = 'true';
else
    nfval = 'false';
end

if s.CenterDC,
    cdval = 'true';
else
    cdval = 'false';
end
s = changedisplay(s, 'NormalizedFrequency', nfval,'CenterDC', cdval);

disp(s);

% [EOF]
