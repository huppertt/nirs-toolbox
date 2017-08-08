function n = nadd(this)
%NADD Returns the number of adders  

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

n = 0;
for k=1:length(this.Stage)
    n = n + nadd(this.Stage(k));
end


% [EOF]
