function verifyautoscalability(this)
%VERIFYAUTOSCALABILITY   

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

for k=1:length(this.Stage)
  verifyautoscalability(this.Stage(k));
end

% [EOF]
