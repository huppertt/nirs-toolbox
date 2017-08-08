function realizeflag = isrealizable(Hd)
%ISREALIZABLE True if the structure can be realized by simulink

%   Author(s): Honglei Chen
%   Copyright 1988-2004 The MathWorks, Inc.

nsections = length(Hd.Stage); 
for k=1:nsections, 
   realizeflag(k) = isrealizable(Hd.Stage(k));
end 
realizeflag = all(realizeflag);
% [EOF]
