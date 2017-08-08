function flag = ishalfnyqinterval(this)
%ISHALFNYQINTERVAL
%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

if strcmpi(this.SpectrumRange,'whole'),
  flag = false;
else
  flag = true;
end


% [EOF]
