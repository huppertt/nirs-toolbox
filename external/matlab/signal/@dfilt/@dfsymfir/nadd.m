function n = nadd(this)
%NADD Returns the number of adders  

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

n = 2*nmult(this,false,false)-1;
if n<0, n = 0; end

% [EOF]
