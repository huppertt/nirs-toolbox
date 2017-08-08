function fdesign = getfdesign(this)
%GETFDESIGN   Get the fdesign.

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

fdesign = this.privfdesign;
if ~isempty(fdesign)
    fdesign = copy(this.privfdesign);
end

% [EOF]
