function flag = isempty(h)
%ISEMPTY   

%   Author(s): V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

flag = isempty(h.Numerator) && isempty(h.Denominator); 


% [EOF]
