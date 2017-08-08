function setfdesign(this, fdesign)
%SETFDESIGN   Set the fdesign.

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

% Copy the FDesign so that users cannot modify what is stored internally.
if ~isempty(fdesign)
    fdesign = copy(fdesign);
end

set(this, 'privfdesign', fdesign, 'privMeasurements', []);

% [EOF]
