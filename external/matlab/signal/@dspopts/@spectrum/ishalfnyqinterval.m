function flag = ishalfnyqinterval(this)
%ISHALFNYQINTERVAL 
%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

if strcmpi(this.SpectrumType,'Twosided'),
  flag = false;
else
  flag = true;
end


% [EOF]
