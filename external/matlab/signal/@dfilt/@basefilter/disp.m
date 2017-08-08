function disp(Hb)
%DISP Object display.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.

if length(this) > 1
    vectordisp(this);
else
    disp(get(Hb))
end


% [EOF]
